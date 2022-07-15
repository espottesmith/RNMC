#ifndef LAT_SOLVER_H
#define LAT_SOLVER_H

#include "../core/sampler.h"
#include <vector>
#include <unordered_map>
#include "../core/solvers.h"
#include <string>
#include <numeric>
#include "LGMC.h"

using namespace LGMC_NS;

struct LatticeUpdate {
    unsigned long int index;
    double propensity;
    int site_one;
    int site_two;
};

struct LatticeEvent {
    int site_one;
    int site_two;
    unsigned long int index;
    double dt;
};

class LatSolver {
    public:
        LatSolver(unsigned long int seed, std::vector<double> &&initial_propensities);
        LatSolver(unsigned long int seed, std::vector<double> &initial_propensities);

        void update(Update update);
        void update(std::vector<Update> updates);

        void update(LatticeUpdate lattice_update);
        void update(std::vector<LatticeUpdate> lattice_updates);

        std::pair<std::optional<Event>, std::optional<LatticeEvent>> event_lattice();

        std::string make_string(int site_one, int site_two);

        double propensity_sum;  

        // for comptability 
        std::optional<Event> event() {return std::optional<Event>();};

    private:
        Sampler sampler;

        int number_of_active_indices;               // end simulation of no sites with non zero propensity            

        unsigned long int last_non_zero_event;   

        std::unordered_map<std::string,                     // lattice propensities 
        std::vector< std::pair<double, int> > > props;      // key: (site_one, site_two) value: propensity 
        std::vector<double> propensities;                   // Gillepsie propensities 

};

/* ---------------------------------------------------------------------- */

// LinearSolver implementation
// LinearSolver can opperate directly on the passed propensities using a move
LatSolver::LatSolver(unsigned long int seed,
    std::vector<double> &&initial_propensities) :
    sampler (Sampler(seed)),
    // if this move isn't here, the semantics is that initial
    // propensities gets moved into a stack variable for the function
    // call and that stack variable is copied into the object.
    propensities (std::move(initial_propensities)),
    number_of_active_indices (0),
    propensity_sum (0.0) 
{
    for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
                last_non_zero_event = i;
            }

        }
};

/* ---------------------------------------------------------------------- */

LatSolver::LatSolver( unsigned long int seed,
    std::vector<double> &initial_propensities) :
    sampler (Sampler(seed)),
    propensities (initial_propensities),
    number_of_active_indices (0),
    propensity_sum (0.0) 
    {
        for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
                last_non_zero_event = i;
            }

        }
    };

/* ---------------------------------------------------------------------- */

void LatSolver::update(Update update) {

    if (propensities[update.index] > 0.0) number_of_active_indices--;

    if (update.propensity > 0.0) {
        number_of_active_indices++;
        if ( update.index > last_non_zero_event )
            last_non_zero_event = update.index;
    }


    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;
    propensities[update.index] = update.propensity;
};

/* ---------------------------------------------------------------------- */

void LatSolver::update(LatticeUpdate lattice_update) {
    
    propensity_sum += lattice_update.propensity;

    std::string hash = make_string(lattice_update.site_one, lattice_update.site_two);
    props[hash].push_back(std::make_pair(lattice_update.propensity, lattice_update.index));
};

/* ---------------------------------------------------------------------- */

void LatSolver::update(std::vector<LatticeUpdate> lattice_updates) {
    for (LatticeUpdate u : lattice_updates) {
        update(u);
    }
};

/* ---------------------------------------------------------------------- */

void LatSolver::update(std::vector<Update> updates) {
    for (Update u : updates) {
        update(u);
    }
};

/* ---------------------------------------------------------------------- */

std::pair<std::optional<Event>, std::optional<LatticeEvent>> LatSolver::event_lattice() {
    if (number_of_active_indices == 0) {
        propensity_sum = 0.0;
        return std::make_pair(std::optional<Event>(), std::optional<LatticeEvent>());
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    unsigned long m;
    bool isFound = false;
    bool isLattice = false;
    unsigned long int reaction_id = 0;
    int site_one = 0;
    int site_two = 0;
    std::string hash;


    // start with Gillespie propensities
    for (m = 0; m < propensities.size(); m++) {
        partial += propensities[m];
        if (partial > fraction) {
            isFound = true;
            reaction_id = m;
            break;
        }
    }

    // go through lattice propensities if not found 
    if(!isFound) {
        auto it = props.begin();
        while(!isFound || it != props.end()) {
            for(int i = 0; i < it->second.size(); i++ ) {
                partial += it->second[i].first;

                if(partial > fraction) {
                    isFound = true;
                    isLattice = true;
                    hash = it->first;
                    reaction_id = it->second[i].second;

                    std::size_t pos = hash.find(".");
                    site_one = stoi(hash.substr(0, pos));
                    site_two = stoi(hash.substr(pos));
                    isFound = true;
                    break;
                }
            }
            it++;
        }
    }
    
    double dt = - std::log(r2) / propensity_sum;

    if(isFound && isLattice) {
        return std::make_pair(std::optional<Event>(), std::optional<LatticeEvent> ( LatticeEvent {.index = reaction_id,
                                                                                .site_one = site_one, .site_two = site_two, 
                                                                                .dt = dt}));
    }
    else if(isFound && !isLattice) {
        return std::make_pair(std::optional<Event>(Event {.index = last_non_zero_event, .dt = dt}), std::optional<LatticeEvent> ());
    }
    else {
        return std::make_pair(std::optional<Event> (Event {.index = last_non_zero_event, .dt = dt}), std::optional<LatticeEvent> ());
    }
        
}

/* ---------------------------------------------------------------------- */

std::string LatSolver::make_string(int site_one, int site_two) {

    return (site_one < site_two) ? std::to_string(site_one) + "." + 
    std::to_string(site_two) : std::to_string(site_two) + "." + std::to_string(site_one);

} // make_string

#endif