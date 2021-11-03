#pragma once
#include <queue>
#include <mutex>
#include <thread>
#include "simulation.h"


struct SeedQueue {
    std::queue<unsigned long int> seeds;
    std::mutex mutex;

    SeedQueue(unsigned long int number_of_seeds, unsigned long int base_seed) {
        for (unsigned long int i = base_seed;
             i < number_of_seeds + base_seed;
             i++) {
            seeds.push(i);
        }
    }

    std::optional<unsigned long int> get_seed() {
        std::lock_guard<std::mutex> lock (mutex);

        if (seeds.empty()) {
            return std::optional<unsigned long int> ();
        } else {
            unsigned long int result = seeds.front();
            seeds.pop();
            return std::optional<unsigned long int> (result);
        }
    }
};

struct HistoryPacket {
    std::vector<HistoryElement> history;
    unsigned long int seed;
};

struct HistoryQueue {
    // the flow of trajectory histories from the simulator threads to
    // the dispatcher is subtle and important. The vector of histories
    // is allocated by the simulator. Once the simulation is finished,
    // it is moved into the history queue. Then the dispatcher moves
    // it out of the history queue, writes it into the initial state
    // database and then frees it. Without all the carefully placed
    // std::move calls, there will be a lot of unnecessary allocations
    // and frees. In typical C++ fashion, this is made significantly
    // more complex by the subtle move semantics of std::optional. If
    // you are going to change this code, you need to use gdb to check
    // that you haven't accidently introduced extra allocations and
    // frees (i.e the vector allocated by the simulation thread points
    // to exactly the same memory as the vector which is used to write
    // to the initial state database, and every time a move is
    // supposed to happen, the old reference is actually zerod out).
    std::queue<HistoryPacket> history_packets;
    std::mutex mutex;

    void insert_history(HistoryPacket history_packet) {
        std::lock_guard<std::mutex> lock (mutex);
        history_packets.push(std::move(history_packet));
    }

    std::optional<HistoryPacket> get_history() {
        std::lock_guard<std::mutex> lock (mutex);
        if (history_packets.empty()) {
            return std::optional<HistoryPacket> ();
        } else {
            HistoryPacket result = std::move(history_packets.front());
            history_packets.pop();
            return std::optional<HistoryPacket> (std::move(result));

        }
    };

};


template <typename Solver>
struct SimulatorPayload {
    ReactionNetwork &reaction_network;
    HistoryQueue &history_queue;
    SeedQueue &seed_queue;
    int step_cutoff;

    SimulatorPayload(
        ReactionNetwork &reaction_network,
        HistoryQueue &history_queue,
        SeedQueue &seed_queue,
        int step_cutoff
        ) :

            reaction_network (reaction_network),
            history_queue (history_queue),
            seed_queue (seed_queue),
            step_cutoff (step_cutoff) {};

    void run_simulator() {

        while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {

            unsigned long int seed = maybe_seed.value();
            Simulation<Solver> simulation (reaction_network, seed, step_cutoff);
            simulation.execute_steps(step_cutoff);
            history_queue.insert_history(
                std::move(
                    HistoryPacket {
                        .history = std::move(simulation.history),
                        .seed = seed
                        }));
        }
    }
};

template <typename Solver>
struct Dispatcher {
    SqlConnection reaction_database;
    SqlConnection initial_state_database;
    ReactionNetwork reaction_network;
    SqlStatement<TrajectoriesSql> trajectories_stmt;
    SqlWriter<TrajectoriesSql> trajectories_writer;
    HistoryQueue history_queue;
    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    int step_cutoff;
    int number_of_simulations;
    int number_of_threads;

    Dispatcher(
        std::string reaction_database_file,
        std::string initial_state_database_file,
        unsigned long int number_of_simulations,
        unsigned long int base_seed,
        int number_of_threads,
        int step_cutoff,
        int dependency_threshold) :
        reaction_database (
            reaction_database_file,
            SQLITE_OPEN_READWRITE),
        initial_state_database (
            initial_state_database_file,
            SQLITE_OPEN_READWRITE),
        reaction_network (
            reaction_database,
            initial_state_database,
            dependency_threshold),
        trajectories_stmt (initial_state_database),
        trajectories_writer (trajectories_stmt),
        history_queue (),
        seed_queue (number_of_simulations, base_seed),

        // don't want to start threads in the constructor.
        threads (),
        step_cutoff (step_cutoff),
        number_of_simulations (number_of_simulations),
        number_of_threads (number_of_threads)
        {};

    void run_dispatcher();
    void record_simulation_history(HistoryPacket history_packet);
};



template <typename Solver>
void Dispatcher<Solver>::run_dispatcher() {

    threads.resize(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        threads[i] = std::thread (
            [](SimulatorPayload<Solver> payload) {payload.run_simulator();},
            SimulatorPayload<Solver> (
                reaction_network,
                history_queue,
                seed_queue,
                step_cutoff)
            );

    }

    int trajectories_written = 0;
    while (trajectories_written < number_of_simulations) {

        std::optional<HistoryPacket>
            maybe_history_packet = history_queue.get_history();

        if (maybe_history_packet) {
            HistoryPacket history_packet = std::move(maybe_history_packet.value());
            record_simulation_history(std::move(history_packet));
            trajectories_written += 1;
        };
    }

    for (int i = 0; i < number_of_threads; i++) threads[i].join();

    initial_state_database.exec(
        "DELETE FROM trajectories WHERE rowid NOT IN"
        "(SELECT MIN(rowid) FROM trajectories GROUP BY seed, step);");

    std::cerr << time_stamp()
              << "removing duplicate trajectories...\n";


};

template <typename Solver>
void Dispatcher<Solver>::record_simulation_history(HistoryPacket history_packet) {
    // TODO: break into subtransactions for super long simulations
    initial_state_database.exec("BEGIN");
    for (unsigned long int i = 0; i < history_packet.history.size(); i++) {
        trajectories_writer.insert(

            TrajectoriesSql {
                .seed = (int) history_packet.seed,
                .step = (int) i,
                .reaction_id = history_packet.history[i].reaction,
                .time = history_packet.history[i].time
            });
    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time_stamp()
              << "wrote trajectory "
              << history_packet.seed
              << " to database\n";
};