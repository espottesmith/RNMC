import math
import numpy as np
import networkx as nx
import os
import pickle
from copy import deepcopy

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import metal_edge_extender


def visualize_molecule_entry(molecule_entry, path):

    atom_colors = {
        'O' : 'red',
        'H' : 'gray',
        'C' : 'black',
        'Li' : 'purple'
        }

    graph = deepcopy(molecule_entry.graph).to_undirected()

    nx.set_node_attributes(graph, "", 'label')
    nx.set_node_attributes(graph, 'filled', 'style')
    nx.set_node_attributes(graph, "circle", 'shape')
    nx.set_node_attributes(graph, "0.2", 'width')
    nx.set_node_attributes(
        graph,
        dict(enumerate([atom_colors[a] for a in molecule_entry.species])),
        'color')


    charge = molecule_entry.charge
    agraph = nx.nx_agraph.to_agraph(graph)
    if charge != 0:
        agraph.add_node(
            'charge',
            label = str(charge),
            fontsize = '25.0',
            shape = 'box',
            color = 'gray',
            style = 'dashed, rounded')

    agraph.layout()
    agraph.draw(path, format='pdf')




class ReactionNetworkSerializationData:
    """
    main attributes
    species_to_index: a dictionary which maps species names to their
    numerical_index which is used to refer to species during simulations.

    index_to_speies: the reverse dictionary to species_to_index

    species_data: maps species indicies to species data.

    index_to_reaction: a list which assigns indices to each reaction.
    Reactants and products are in terms of species indices.

    dependency_graph: a dictionary mapping each reaction index to
    all the reaction indices which have a reactant in the reactants or products
    of the reaction. In other words, it maps a reaction index to the list
    of all reactions whose propensities need to be updated when the reaction
    occours
    """

    def pp_reaction(self,index):
        """
        pretty print a reaction given its index
        """
        reaction = self.index_to_reaction[index]
        reactants = " + ".join([str(self.index_to_species[reactant_index])
                                for reactant_index in reaction['reactants']])
        products = " + ".join([str(self.index_to_species[product_index])
                                for product_index in reaction['products']])
        dG = str(reaction['free_energy'])
        return (reactants + " -> " + products).ljust(50) + dG

    def visualize_molecules(self):
        folder = self.network_folder + '/molecule_diagrams'
        if os.path.isdir(folder):
            return

        os.mkdir(folder)
        for index in range(self.number_of_species):
            molecule_entry = self.species_data[index]
            visualize_molecule_entry(
                molecule_entry,
            folder + '/' + str(index) + '.pdf')


    def find_index_from_mol_graph(self, mol_graph_file_path, charge):
        """
        given a file 'molecule.xyz', find the index corresponding to the
        molecule graph with given charge
        """
        target_mol_graph = MoleculeGraph.with_local_env_strategy(
            Molecule.from_file(mol_graph_file_path),
            OpenBabelNN())

        ### correction to the molecule graph
        target_mol_graph = metal_edge_extender(target_mol_graph)

        match = False
        species_index = -1
        while not match:
            species_index += 1
            data = self.species_data[species_index]
            species_mol_graph = data.mol_graph

            if data.charge == charge:
                match = target_mol_graph.isomorphic_to(species_mol_graph)

        if match:
            return species_index
        else:
            return None


    def __extract_index_mappings(self,reactions):
        """
        assign each species an index and construct
        forward and backward mappings between indicies and species.

        assign each reaction an index and construct
        a mapping from reaction indices to reaction data
        """
        species_to_index = {}
        index_to_reaction = []
        index = 0
        reaction_count = 0

        for reaction in reactions:
            reaction_count += 1
            entry_ids = {e.entry_id for e in reaction.reactants + reaction.products}
            for entry_id in entry_ids:
                species = entry_id
                if species not in species_to_index:
                    species_to_index[species] = index
                    index = index + 1


            reactant_indices = [species_to_index[reactant]
                                for reactant in reaction.reactant_ids]
            product_indices = [species_to_index[product]
                               for product in reaction.product_ids]

            forward_free_energy = reaction.free_energy_A
            backward_free_energy = reaction.free_energy_B


            index_to_reaction.append({'reactants' : reactant_indices,
                                      'products' : product_indices,
                                      'free_energy' : forward_free_energy})
            index_to_reaction.append({'reactants' : product_indices,
                                      'products' : reactant_indices,
                                      'free_energy' : backward_free_energy})


        for reaction in index_to_reaction:
            dG = reaction['free_energy']
            if dG > 0:
                rate = math.exp(- self.positive_weight_coef * dG)
            else:
                rate = math.exp(- dG)
            reaction['rate'] = rate

        rev = {i : species for species, i in species_to_index.items()}
        self.number_of_reactions = 2 * reaction_count
        self.number_of_species = index
        self.species_to_index = species_to_index
        self.index_to_species = rev
        self.index_to_reaction = index_to_reaction

    def __extract_species_data(self,entries_list):
        species_data = {}
        for entry in entries_list:
            entry_id = entry.entry_id
            if entry_id in self.species_to_index:
                species_data[self.species_to_index[entry_id]] = entry

        self.species_data = species_data

    def __init__(self,reaction_network,
                 # list of triples (path to molecule.xyz, charge, number)
                 initial_state_data,
                 network_folder,
                 param_folder,
                 logging = False,
                 positive_weight_coef = 39):
        """
        Input: a reaction network object
        """
        reactions = reaction_network
        entries_list = reaction_network.entries_list
        self.network_folder = network_folder
        self.param_folder = param_folder
        self.logging = logging
        self.positive_weight_coef = positive_weight_coef


        self.__extract_index_mappings(reactions)
        if logging:
            print("extracted index mappings")

        self.__extract_species_data(entries_list)
        if logging:
            print("extracted species data")

        self.initial_state = np.zeros(self.number_of_species)
        for (path, charge, count) in initial_state_data:
            index = self.find_index_from_mol_graph(path, charge)
            self.initial_state[index] = count

        if logging:
            print("set initial state")

        if logging:
            print("finished building serialization data")




def serialize_reaction_network(rnsd):
    number_of_species_postfix = "/number_of_species"
    number_of_reactions_postfix = "/number_of_reactions"
    number_of_reactants_postfix = "/number_of_reactants"
    reactants_postfix = "/reactants"
    number_of_products_postfix = "/number_of_products"
    products_postfix = "/products"
    factor_zero_postfix = "/factor_zero"
    factor_two_postfix = "/factor_two"
    factor_duplicate_postfix = "/factor_duplicate"
    rates_postfix = "/rates"
    initial_state_postfix = "/initial_state"

    folder = rnsd.network_folder

    os.mkdir(folder)

    with open(folder + number_of_species_postfix, 'w') as f:
        f.write(str(rnsd.number_of_species) + '\n')

    with open(folder + number_of_reactions_postfix, 'w') as f:
        f.write(str(rnsd.number_of_reactions) + '\n')

    with open(folder + number_of_reactants_postfix, 'w') as f:
        for reaction in rnsd.index_to_reaction:
            f.write(str(len(reaction['reactants'])) + '\n')

    with open(folder + reactants_postfix, 'w') as f:
        for reaction in rnsd.index_to_reaction:
            for index in reaction['reactants']:
                f.write(str(index) + ' ')
            f.write('\n')

    with open(folder + number_of_products_postfix, 'w') as f:
        for reaction in rnsd.index_to_reaction:
            f.write(str(len(reaction['products'])) + '\n')

    with open(folder + products_postfix, 'w') as f:
        for reaction in rnsd.index_to_reaction:
            for index in reaction['products']:
                f.write(str(index) + ' ')
            f.write('\n')

    with open(folder + factor_two_postfix, 'w') as f:
        f.write(('%e' % 1.0) + '\n')

    with open(folder + factor_zero_postfix, 'w') as f:
        f.write(('%e' % 1.0) + '\n')

    with open(folder + factor_duplicate_postfix, 'w') as f:
        f.write(('%e' % 1.0) + '\n')

    with open(folder + rates_postfix, 'w') as f:
        for reaction in rnsd.index_to_reaction:
            f.write(('%e' % reaction['rate']) + '\n')

    with open(folder + initial_state_postfix, 'w') as f:
        for i in range(rnsd.number_of_species):
            f.write(str(int(rnsd.initial_state[i])) + '\n')

    with open(folder + "/rnsd.pickle",'wb') as f:
        pickle.dump(rnsd, f)

    print("finished serializing")


def serialize_simulation_parameters(seeds,
                                    number_of_threads,
                                    time_cutoff,
                                    folder):

    number_of_seeds_postfix = "/number_of_seeds"
    number_of_threads_postfix = "/number_of_threads"
    seeds_postfix = "/seeds"
    time_cutoff_postfix = "/time_cutoff"

    os.mkdir(folder)

    with open(folder + number_of_seeds_postfix, 'w') as f:
        f.write(str(len(seeds)) + '\n')

    with open(folder + number_of_threads_postfix, 'w') as f:
        f.write(str(number_of_threads) + '\n')

    with open(folder + seeds_postfix, 'w') as f:
        for seed in seeds:
            f.write(str(seed) + '\n')

    with open(folder + time_cutoff_postfix, 'w') as f:
        f.write(('%e' % time_cutoff) + '\n')


def collect_duplicate_pathways(pathways):
    pathway_dict = {}
    for pathway in pathways:
        key = frozenset(pathway)
        if key in pathway_dict:
            pathway_dict[key]['frequency'] += 1
        else:
            pathway_dict[key] = {'pathway' : pathway, 'frequency' : 1}
    return pathway_dict

def update_state(state, reaction):
    for species_index in reaction['reactants']:
        state[species_index] -= 1

    for species_index in reaction['products']:
        state[species_index] += 1

class SimulationAnalyser:
    def extract_reaction_pathways(self, target_species_index):
        """
        given a reaction history and a target molecule, find the
        first reaction which produced the target molecule (if any).
        Apply that reaction to the initial state to produce a partial
        state array. Missing reactants have negative values in the
        partial state array. Now loop through the reaction history
        to resolve the missing reactants.
        """
        reaction_pathway_list = []
        for reaction_history in self.reaction_histories:

            # -1 if target wasn't produced
            # index of reaction if target was produced
            reaction_producing_target_index = -1
            for reaction_index in reaction_history:
                reaction = self.rnsd.index_to_reaction[reaction_index]
                if target_species_index in reaction['products']:
                    reaction_producing_target_index = reaction_index
                    break

            if reaction_producing_target_index == -1:
                continue
            else:
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


    def pp_pathways(self,target_species_index):
        if target_species_index not in self.reaction_pathways_dict:
            self.extract_reaction_pathways(target_species_index)

        pathways = self.reaction_pathways_dict[target_species_index]

        for _, unique_pathway in sorted(
                pathways.items(),
                key = lambda item: -item[1]['frequency']):

            print(str(unique_pathway['frequency']) + " occurrences:")
            for reaction_index in unique_pathway['pathway']:
                print(self.rnsd.pp_reaction(reaction_index))
            print()

    def generate_pathway_report(self, target_species_index):
        self.rnsd.visualize_molecules()
        folder = self.rnsd.network_folder + '/report_' + str(target_species_index)
        os.mkdir(folder)

        with open(folder + '/report.tex','w') as f:
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



    def __init__(self, rnsd, network_folder):

        history_postfix = '/simulation_histories'
        simulation_histories_folder = network_folder + history_postfix
        self.rnsd = rnsd
        self.initial_state = rnsd.initial_state
        self.reaction_pathways_dict = {}
        self.reaction_histories = []

        for filename in os.listdir(simulation_histories_folder):
            history = []
            with open(simulation_histories_folder + '/' + filename, 'r') as f:
                for line in f:
                    history.append(int(line))

            self.reaction_histories.append(history)
