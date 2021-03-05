import numpy as np
import os
import pickle
import copy
import math
from typing import Tuple, Optional, Union, List, Dict
from pathlib import Path
import random
import sys

from pymatgen.analysis.graphs import MoleculeGraph

from mrnet.core.reactions import Reaction
from mrnet.network.reaction_network import ReactionNetwork
from mrnet.network.reaction_generation import ReactionGenerator

from rnmc.visualize import visualize_molecule_entry


class SerializedReactionNetwork:
    """
    An object designed to store data from a ReactionNetwork suitable for use with
    the C RNMC code.
    """

    def __init__(self,
                 reaction_network: Union[ReactionNetwork, ReactionGenerator],
                 initial_state_data: List[Tuple[MoleculeGraph, Union[int, float]]],
                 network_folder: Path,
                 param_folder: Path,
                 logging: bool = False,
                 positive_weight_coefficient: float = 38.61):
        """
        Args:
            reaction_network (ReactionNetwork):
            initial_state_data (tuple (MoleculeGraph, count)):
            network_folder (Path):
            param_folder (Path):
            logging (bool):
            positive_weight_coefficient (float):
        """
        if isinstance(reaction_network, ReactionGenerator):
            reactions = reaction_network
        else:
            reactions = reaction_network.reactions
        entries_list = reaction_network.entries_list
        self.network_folder = network_folder
        self.param_folder = param_folder
        self.logging = logging
        self.positive_weight_coefficient = positive_weight_coefficient

        self._extract_index_mappings(reactions)
        if self.logging:
            print("extracted index mappings")

        self._extract_species_data(entries_list)
        if self.logging:
            print("extracted species data")

        self.initial_state = np.zeros(self.number_of_species)
        for (mg, count) in initial_state_data:
            index = self.find_index_from_mol_graph(mg)
            self.initial_state[index] = count

        if self.logging:
            print("set initial state")

        if self.logging:
            print("finished building serialization data")

    def pp_reaction(self,
                    index: int
                    ) -> str:
        """
        pretty print a reaction given its index

        Args:
            index (int): Index of the reaction to be printed

        Return:
            str (formatted, pretty-printed reaction representation)
        """
        reaction = self.index_to_reaction[index]
        reactants = " + ".join([str(self.index_to_species[reactant_index])
                                for reactant_index in reaction['reactants']])
        products = " + ".join([str(self.index_to_species[product_index])
                                for product_index in reaction['products']])
        dG = str(reaction['free_energy'])
        return (reactants + " -> " + products).ljust(50) + dG

    def visualize_molecules(self):
        """
        Visualize molecules in the serialized reaction network

        Args:
            None

        Return:
            None
        """
        folder = self.network_folder / 'molecule_diagrams'
        if folder.is_dir():
            return

        folder.mkdir(exist_ok=True)
        for index in range(self.number_of_species):
            molecule_entry = self.species_data[index]
            visualize_molecule_entry(
                molecule_entry,
                folder / (str(index) + '.pdf')
            )

    def find_index_from_mol_graph(self,
                                  target_mol_graph: MoleculeGraph) -> Optional[int]:
        """
        Given MoleculeGraph, find the index of the corresponding entry in the
        SerializedReactionNetwork species_data.

        Args:
            target_mol_graph: MoleculeGraph of interest

        Return:
            index: index of the molecule matched by target_mol_graph
                (or None, if no match is found)
        """

        match = False
        for index, data in self.species_data.items():
            species_mol_graph = data.mol_graph
            if data.charge == target_mol_graph.molecule.charge:
                match = target_mol_graph.isomorphic_to(species_mol_graph)
            if match:
                return index

        return None

    def _extract_index_mappings(self,
                                reactions: List[Reaction],
                                use_thermo_cost: bool = False):
        """
        Assign each species an index and construct
        forward and backward mappings between indicies and species.

        Args:
            reactions (list of Reaction objects): reactions to be used for mapping

        Return:
            None

        """
        species_to_index = dict()
        index_to_reaction = list()
        index = 0
        reaction_count = 0

        for reaction in reactions:
            reaction_count += 1
            entry_ids = {e.entry_id for e in reaction.reactants + reaction.products}
            for entry_id in entry_ids:
                species = entry_id
                if species not in species_to_index:
                    species_to_index[species] = index
                    index += 1

            reactant_indices = [species_to_index[reactant]
                                for reactant in reaction.reactant_ids]
            product_indices = [species_to_index[product]
                               for product in reaction.product_ids]

            forward_free_energy = reaction.free_energy_A
            backward_free_energy = reaction.free_energy_B

            forward_rate = reaction.k_A
            backward_rate = reaction.k_B

            index_to_reaction.append({'reactants' : reactant_indices,
                                      'products' : product_indices,
                                      'free_energy' : forward_free_energy,
                                      'rate_constant': forward_rate})
            index_to_reaction.append({'reactants' : product_indices,
                                      'products' : reactant_indices,
                                      'free_energy' : backward_free_energy,
                                      'rate_constant': backward_rate})

        if use_thermo_cost:
            for reaction in index_to_reaction:
                dG = reaction['free_energy']
                if dG > 0:
                    rate = math.exp(- self.positive_weight_coefficient * dG)
                else:
                    rate = math.exp(- dG)
                reaction['rate_constant'] = rate

        rev = {i : species for species, i in species_to_index.items()}
        self.number_of_reactions = 2 * reaction_count
        self.number_of_species = index
        self.species_to_index = species_to_index
        self.index_to_species = rev
        self.index_to_reaction = index_to_reaction

    def _extract_species_data(self, entries_list):
        species_data = dict()
        for entry in entries_list:
            entry_id = entry.entry_id
            if entry_id in self.species_to_index:
                species_data[self.species_to_index[entry_id]] = entry

        self.species_data = species_data

    def serialize(self):
        """
        Dump data into files

        Args:
            None
        """
        folder = self.network_folder

        folder.mkdir(exist_ok=True)

        with open((folder / "number_of_species").as_posix(), 'w') as f:
            f.write(str(self.number_of_species) + '\n')

        with open((folder / "number_of_reactions").as_posix(), 'w') as f:
            f.write(str(self.number_of_reactions) + '\n')

        with open((folder / "number_of_reactants").as_posix(), 'w') as f:
            for reaction in self.index_to_reaction:
                f.write(str(len(reaction['reactants'])) + '\n')

        with open((folder / "reactants").as_posix(), 'w') as f:
            for reaction in self.index_to_reaction:
                for index in reaction['reactants']:
                    f.write(str(index) + ' ')
                f.write('\n')

        with open((folder / "number_of_products").as_posix(), 'w') as f:
            for reaction in self.index_to_reaction:
                f.write(str(len(reaction['products'])) + '\n')

        with open((folder / "products").as_posix(), 'w') as f:
            for reaction in self.index_to_reaction:
                for index in reaction['products']:
                    f.write(str(index) + ' ')
                f.write('\n')

        with open((folder / "factor_two").as_posix(), 'w') as f:
            f.write(('%e' % 1.0) + '\n')

        with open((folder / "factor_zero").as_posix(), 'w') as f:
            f.write(('%e' % 1.0) + '\n')

        with open((folder / "factor_duplicate").as_posix(), 'w') as f:
            f.write(('%e' % 1.0) + '\n')

        with open((folder / "rates").as_posix(), 'w') as f:
            for reaction in self.index_to_reaction:
                f.write(('%e' % reaction['rate_constant']) + '\n')

        with open((folder / "initial_state").as_posix(), 'w') as f:
            for i in range(self.number_of_species):
                f.write(str(int(self.initial_state[i])) + '\n')

        with open((folder / "rnsd.pickle").as_posix(),'wb') as f:
            pickle.dump(self, f)

        print("finished serializing")


class SimulationAnalyser:
    """
    A class to analyze the resutls of a set of MC runs
    """

    def __init__(self,
                 rnsd: SerializedReactionNetwork,
                 network_folder: Path):
        """
        Params:
            rnsd (SerializedReactionNetwork):
            network_folder (Path):
        """

        self.network_folder = network_folder
        self.histories_folder =  network_folder / 'simulation_histories'
        self.rnsd = rnsd
        self.initial_state = rnsd.initial_state
        self.reaction_pathways_dict = dict()
        self.reaction_histories = list()
        self.time_histories = list()

        histories_contents = sorted([x.name for x in self.histories_folder.iterdir()])
        reaction_histories_contents = [x for x in histories_contents if x.startswith("reactions")]
        time_histories_contents = [x for x in histories_contents if x.startswith("times")]

        reaction_seeds = [x.split("_")[1] for x in reaction_histories_contents]
        time_seeds = [x.split("_")[1] for x in reaction_histories_contents]

        if reaction_seeds != time_seeds:
            raise ValueError("Reactions and times not from same set of initial seeds!")

        for filename in reaction_histories_contents:
            history = list()
            with open((self.histories_folder / filename).as_posix()) as f:
                for line in f:
                    history.append(int(line.strip()))

            self.reaction_histories.append(np.array(history))

        for filename in time_histories_contents:
            history = list()
            with open((self.histories_folder / filename).as_posix()) as f:
                for line in f:
                    history.append(float(line.strip()))

            self.time_histories.append(np.array(history))

        self.number_simulations = len(self.reaction_histories)

    def extract_reaction_pathways(self,
                                  target_species_index: int):
        """
        given a reaction history and a target molecule, find the
        first reaction which produced the target molecule (if any).
        Apply that reaction to the initial state to produce a partial
        state array. Missing reactants have negative values in the
        partial state array. Now loop through the reaction history
        to resolve the missing reactants.

        Args:
            target_species_index (int): Product molecule of interest
        """
        reaction_pathway_list = list()
        for reaction_history in self.reaction_histories:

            reaction_producing_target_index = None
            for reaction_index in reaction_history:
                reaction = self.rnsd.index_to_reaction[reaction_index]
                if target_species_index in reaction['products']:
                    reaction_producing_target_index = reaction_index
                    break

            if reaction_producing_target_index is None:
                continue

            pathway = [reaction_producing_target_index]
            partial_state = np.copy(self.initial_state)
            final_reaction = self.rnsd.index_to_reaction[pathway[0]]
            update_state(partial_state, final_reaction)

            negative_species = list(np.where(partial_state < 0)[0])

            while(len(negative_species) != 0):
                for species_index in negative_species:
                    for reaction_index in reaction_history:
                        reaction = self.rnsd.index_to_reaction[reaction_index]
                        if species_index in reaction['products']:
                            update_state(partial_state, reaction)
                            pathway.insert(0,reaction_index)
                            break

                negative_species = list(np.where(partial_state < 0)[0])

            reaction_pathway_list.append(pathway)

        reaction_pathway_dict = collect_duplicate_pathways(reaction_pathway_list)
        self.reaction_pathways_dict[target_species_index] = reaction_pathway_dict

    def pp_pathways(self,
                    target_species_index: int):
        """
        Pretty-print pathways to a molecule of interest

        Args:
            target_species_index (int): Index of the product molecule of interest

        Returns:
            None
        """

        if target_species_index not in self.reaction_pathways_dict:
            self.extract_reaction_pathways(target_species_index)

        pathways = self.reaction_pathways_dict[target_species_index]

        for _, unique_pathway in sorted(
                pathways.items(),
                key=lambda item: -item[1]['frequency']):

            print(str(unique_pathway['frequency']) + " occurrences:")
            for reaction_index in unique_pathway['pathway']:
                print(self.rnsd.pp_reaction(reaction_index))
            print()

    def generate_pathway_report(self,
                                target_species_index: int):
        """
        Generate a LaTeX document containing pathways to a molecule of interest

        Args:
            target_species_index (int): Index of the product molecule of interest

        Returns:
            None
        """
        self.rnsd.visualize_molecules()
        folder = self.network_folder / ('report_' + str(target_species_index))
        folder.mkdir(exist_ok=True)

        with open((folder / 'report.tex').as_posix(),'w') as f:
            if target_species_index not in self.reaction_pathways_dict:
                self.extract_reaction_pathways(target_species_index)

            pathways = self.reaction_pathways_dict[target_species_index]

            f.write('\\documentclass{article}\n')
            f.write('\\usepackage{graphicx}')
            f.write('\\usepackage[margin=1cm]{geometry}')
            f.write('\\usepackage{amsmath}')
            f.write('\\pagenumbering{gobble}')
            f.write('\\begin{document}\n')

            for _, unique_pathway in sorted(
                    pathways.items(),
                    key = lambda item: -item[1]['frequency']):

                f.write(str(unique_pathway['frequency']) +
                        " occurrences:\n")

                for reaction_index in unique_pathway['pathway']:
                    reaction = self.rnsd.index_to_reaction[reaction_index]
                    f.write('$$\n')
                    first = True
                    for reactant_index in reaction['reactants']:
                        if first:
                            first = False
                        else:
                            f.write('+\n')

                        f.write(
                            '\\raisebox{-.5\\height}{'
                            + '\\includegraphics[scale=0.2]{../molecule_diagrams/'
                            + str(reactant_index)
                            + '.pdf}}\n')

                    f.write('\\xrightarrow{'
                            + ('%.2f' % reaction['free_energy'])
                            + '}\n')

                    first = True
                    for product_index in reaction['products']:
                        if first:
                            first = False
                        else:
                            f.write('+\n')

                        f.write(
                            '\\raisebox{-.5\\height}{'
                            + '\\includegraphics[scale=0.2]{../molecule_diagrams/'
                            + str(product_index)
                            + '.pdf}}\n')

                    f.write('$$')
                    f.write('\n\n\n')

                f.write('\\newpage\n')

            f.write('\\end{document}')

    def generate_time_dependent_profiles(self) -> Dict:
        """
        Generate plottable time-dependent profiles of species and rxns from raw KMC output, obtain final states.

        Args:
            None

        Returns:
            dict containing species profiles, reaction profiles, and final states from each simulation.
                {species_profiles: [ {mol_ind1: [(t0, n(t0)), (t1, n(t1)...], mol_ind2: [...] ,  ... }, {...}, ... ]
                reaction_profiles: [ {rxn_ind1: [t0, t1, ...], rxn_ind2: ..., ...}, {...}, ...]
                final_states: [ {mol_ind1: n1, mol_ind2: ..., ...}, {...}, ...] }

        """
        species_profiles = list()
        reaction_profiles = list()
        final_states = list()

        for n_sim in range(self.number_simulations):
            cumulative_time = self.time_histories[n_sim]
            sim_rxn_history = self.reaction_histories[n_sim]
            sim_species_profile = dict()
            sim_rxn_profile = dict()
            state = dict()
            for index, initial_value in enumerate(self.rnsd.initial_state):
                sim_species_profile[index] = [(0.0, initial_value)]
                state[index] = initial_value
            total_iterations = len(sim_rxn_history)

            for iter in range(total_iterations):
                rxn_ind = sim_rxn_history[iter]
                t = cumulative_time[iter]
                if rxn_ind not in sim_rxn_profile:
                    sim_rxn_profile[rxn_ind] = [t]
                else:
                    sim_rxn_profile[rxn_ind].append(t)

                reacts = self.rnsd.index_to_reaction[rxn_ind]["reactants"]
                prods = self.rnsd.index_to_reaction[rxn_ind]["products"]

                for r_ind in reacts:
                    if r_ind == -1:
                        continue
                    try:
                        state[r_ind] -= 1
                        if state[r_ind] < 0:
                            raise ValueError(
                                "State invalid: negative specie: {}".format(r_ind)
                            )
                        sim_species_profile[r_ind].append((t, state[r_ind]))
                    except KeyError:
                        raise ValueError(
                            "Reactant specie {} given is not in state!".format(
                                r_ind
                            )
                        )
                for p_ind in prods:
                    if p_ind == -1:
                        continue
                    else:
                        if (p_ind in state) and (p_ind in sim_species_profile):
                            state[p_ind] += 1
                            sim_species_profile[p_ind].append((t, state[p_ind]))
                        else:
                            state[p_ind] = 1
                            sim_species_profile[p_ind] = [(0.0, 0), (t, state[p_ind])]

            # for plotting convenience, add data point at final time
            for mol_ind in sim_species_profile:
                sim_species_profile[mol_ind].append(
                    (cumulative_time[-1], state[mol_ind])
                )

            species_profiles.append(sim_species_profile)
            reaction_profiles.append(sim_rxn_profile)
            final_states.append(state)

        return {
            "species_profiles": species_profiles,
            "reaction_profiles": reaction_profiles,
            "final_states": final_states,
        }

    def generate_final_states(self) -> List:
        """
        Generate plottable time-dependent profiles of species and rxns from raw KMC output, obtain final states.

        Args:
            None

        Returns:
            dict containing species profiles, reaction profiles, and final states from each simulation.
                {species_profiles: [ {mol_ind1: [(t0, n(t0)), (t1, n(t1)...], mol_ind2: [...] ,  ... }, {...}, ... ]
                reaction_profiles: [ {rxn_ind1: [t0, t1, ...], rxn_ind2: ..., ...}, {...}, ...]
                final_states: [ {mol_ind1: n1, mol_ind2: ..., ...}, {...}, ...] }

        """
        final_states = list()

        for n_sim in range(self.number_simulations):
            sim_rxn_history = self.reaction_histories[n_sim]
            state = dict()
            for index, initial_value in enumerate(self.rnsd.initial_state):
                state[index] = initial_value
            total_iterations = len(sim_rxn_history)

            for iter in range(total_iterations):
                rxn_ind = sim_rxn_history[iter]

                reacts = self.rnsd.index_to_reaction[rxn_ind]["reactants"]
                prods = self.rnsd.index_to_reaction[rxn_ind]["products"]

                for r_ind in reacts:
                    if r_ind == -1:
                        continue
                    try:
                        state[r_ind] -= 1
                        if state[r_ind] < 0:
                            raise ValueError(
                                "State invalid: negative specie: {}".format(r_ind)
                            )
                    except KeyError:
                        raise ValueError(
                            "Reactant specie {} given is not in state!".format(
                                r_ind
                            )
                        )
                for p_ind in prods:
                    if p_ind == -1:
                        continue
                    else:
                        if p_ind in state:
                            state[p_ind] += 1
                        else:
                            state[p_ind] = 1

            final_states.append(state)

        return final_states

    def final_state_analysis(self, final_states: List[Dict[int, Union[int, float]]]) -> List[Tuple]:
        """
        Gather statistical analysis of the final states of simulation.

        Args:
            final_states: list of dicts of final states, as generated in generate_time_dep_profiles()

        Return:
            sorted_analyzed_states: list of tuples containing statistical data for each species,
                sorted from highest to low avg occurrence
        """
        state_arrays = dict()  # For each molecule, compile an array of its final amounts
        for iter, final_state in enumerate(final_states):
            for mol_ind, amt in final_state.items():
                # Store the amount, and convert key from mol_ind to entry_id
                if self.rnsd.index_to_species[mol_ind] not in state_arrays:
                    state_arrays[self.rnsd.index_to_species[mol_ind]] = np.zeros(
                        self.number_simulations
                    )
                state_arrays[self.rnsd.index_to_species[mol_ind]][iter] = amt
        analyzed_states = dict()  # will contain statistical results of final states
        for mol_entry, state_array in state_arrays.items():
            analyzed_states[mol_entry] = (np.mean(state_array), np.std(state_array))
        # Sort from highest avg final amount to lowest
        sorted_analyzed_states = sorted(
            [(entry_id, data_tup) for entry_id, data_tup in analyzed_states.items()],
            key=lambda x: x[1][0],
            reverse=True,
        )
        return sorted_analyzed_states

    def quantify_rank_reactions(self) -> List[Tuple]:
        """
        Given reaction histories, identify the most commonly occurring reactions, on average.
        Can rank generally, or by reactions of a certain type.

        Args:
            reaction_profiles (list of dicts): reactions fired as a function of time
            reaction_type (string)
            num_rxns (int): the amount of reactions interested in collecting data on. If None, record for all.

        Returns:
            reaction_data: list of reactions and their avg, std of times fired. Sorted by the average times fired.
            [(rxn1, (avg, std)), (rxn2, (avg, std)) ... ]
        """

        reaction_data = dict()  # keeping record of each iteration
        # Loop to count all reactions fired
        for n_sim in range(self.number_simulations):
            rxns_fired = set(self.reaction_histories[n_sim])
            relevant_rxns = rxns_fired

            for rxn_ind in relevant_rxns:
                if rxn_ind not in reaction_data:
                    reaction_data[rxn_ind] = list()
                reaction_data[rxn_ind].append(
                    np.sum(self.reaction_histories[n_sim] == rxn_ind)
                )

        reaction_analysis = dict()
        for rxn_ind, counts in reaction_data.items():
            reaction_analysis[rxn_ind] = (
                np.mean(np.array(counts)),
                np.std(np.array(counts)),
            )

        # Sort reactions by the average amount fired
        sorted_reaction_analysis = sorted(
            [(i, c) for i, c in reaction_analysis.items()],
            key=lambda x: x[1][0],
            reverse=True,
        )
        return sorted_reaction_analysis


def serialize_simulation_parameters(
        folder: Path,
        number_of_threads: int,
        step_cutoff: Optional[int] = 1000000,
        time_cutoff: Optional[float] = None,
        seeds: Optional[List[int]] = None,
        number_of_simulations: Optional[int] = 100,
):
    """

    Args:
        folder (Path): Folder in which to store simulation parameters
        number_of_threads (int): Number of threads to use in simulation
        step_cutoff (int, or None): Number of steps to allow in each simulation.
            Default is 1,000,000
        time_cutoff (float, or None): Time duration of each simulation. Default
            is None
        seeds (List of ints, or None): Random seeds to use for simulations.
            Default is None, meaning that these seeds will be randomly
            generated.
        number_of_simulations (int, or None): Number of simulations. If seeds is
            None (default), then this number (default 100) of random seeds will
            be randomly generated.
    """

    if seeds is not None:
        number_of_seeds = len(seeds)
        random_seeds = seeds
    elif number_of_simulations is not None:
        number_of_seeds = number_of_simulations
        random_seeds = random.sample(list(range(1, sys.maxsize)), number_of_simulations)
    else:
        raise ValueError("Need either number of simulations or set of seeds to proceed!")

    folder.mkdir(exist_ok=True)

    if step_cutoff is not None:
        with open((folder / "step_cutoff").as_posix(), 'w') as f:
            f.write(('%d' % step_cutoff) + '\n')
    elif time_cutoff is not None:
        with open((folder / "time_cutoff").as_posix(), 'w') as f:
            f.write(('%f' % time_cutoff) + '\n')
    else:
        raise ValueError("Either time_cutoff or step_cutoff must be set!")

    folder.mkdir(exist_ok=True)

    with open((folder / "number_of_seeds").as_posix(), 'w') as f:
        f.write(str(number_of_seeds) + '\n')

    with open((folder / "number_of_threads").as_posix(), 'w') as f:
        f.write(str(number_of_threads) + '\n')

    with open((folder / "seeds").as_posix(), 'w') as f:
        for seed in random_seeds:
            f.write(str(seed) + '\n')


def collect_duplicate_pathways(pathways: List[List[int]]) -> Dict:
    pathway_dict = {}
    for pathway in pathways:
        key = frozenset(pathway)
        if key in pathway_dict:
            pathway_dict[key]['frequency'] += 1
        else:
            pathway_dict[key] = {'pathway' : pathway, 'frequency' : 1}
    return pathway_dict


def update_state(
        state: Dict[int, Union[int, float]],
        reaction: Dict):
    for species_index in reaction['reactants']:
        state[species_index] -= 1

    for species_index in reaction['products']:
        state[species_index] += 1