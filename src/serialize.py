import sys
sys.path.append('./mrnet/src')
from mrnet.network.reaction_network import ReactionNetwork
from monty.serialization import loadfn
from rnmc import *


if len(sys.argv) != 4:
    print("usage: python serialize.py json network_folder param_folder")
    exit()

molecule_list_json = sys.argv[1]
network_folder = sys.argv[2]
param_folder = sys.argv[3]

molecule_entries = loadfn(molecule_list_json)
reaction_network = ReactionNetwork.from_input_entries(molecule_entries)
reaction_network.build()

foo = ReactionNetwork.identify_concerted_rxns_via_intermediates(
    reaction_network,
    [e.parameters["ind"] for e in reaction_network.entries_list])

concerted_reactions = ReactionNetwork.add_concerted_rxns(
    reaction_network,
    reaction_network,
    foo)

reaction_network.reactions.extend(concerted_reactions)

initial_state_data = [
    ('./Li.xyz', 1, 30),
    ('./mrnet/test_files/reaction_network_files/EC.xyz',0,30)
    ]

rnsd = ReactionNetworkSerializationData(reaction_network,
                                        initial_state_data,
                                        network_folder,
                                        param_folder,
                                        logging = True)



serialize_reaction_network(rnsd)
serialize_simulation_parameters(seeds=range(1000,100000),
                                number_of_threads=7,
                                time_cutoff=5.0,
                                folder=rnsd.param_folder)

print("run simulation with ../RNMC " + network_folder + " " + param_folder )
