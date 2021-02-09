import sys
sys.path.append('./mrnet/src')
from mrnet.network.reaction_generation import *
from monty.serialization import loadfn
from RNMC import *


if len(sys.argv) != 4:
    print("usage: python serialize.py json network_folder param_folder")
    exit()

molecule_list_json = sys.argv[1]
network_folder = sys.argv[2]
param_folder = sys.argv[3]

molecule_entries = loadfn(molecule_list_json)
reaction_generator = ReactionGenerator(molecule_entries)


initial_state_data = [
    ('./Li.xyz', 1, 30),
    ('./mrnet/test_files/reaction_network_files/EC.xyz',0,30)
    ]

rnsd = ReactionNetworkSerializationData(reaction_generator,
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
