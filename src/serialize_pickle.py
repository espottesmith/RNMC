import sys
sys.path.append('./mrnet/src')
import pickle
import random

from mrnet.network.reaction_network import *
from RNMC import *


if len(sys.argv) != 4:
    print("usage: python serialize.py json network_folder param_folder")
    exit()

network_pickle = sys.argv[1]
network_folder = sys.argv[2]
param_folder = sys.argv[3]

reaction_network = None
with open(network_pickle, "rb") as pick:
    reaction_network = pickle.load(pick)

initial_state_data = [
    ('./Li.xyz', 1, 600),
    ('./mrnet/test_files/reaction_network_files/EC.xyz', 0, 7850),
    ('./mrnet/test_files/reaction_network_files/H2O.xyz', 0, 10),
    ('./mrnet/test_files/reaction_network_files/H.xyz', 1, 10),
    ('./mrnet/test_files/reaction_network_files/OH.xyz', -1, 10),
    ('./CO2.xyz', 0, 50)
    ]

rnsd = ReactionNetworkSerializationData(reaction_network,
                                        initial_state_data,
                                        network_folder,
                                        param_folder,
                                        logging = True)
serialize_reaction_network(rnsd)

if not os.path.exists(rnsd.param_folder):
    serialize_simulation_parameters(seeds=random.sample(list(range(1, 1000000000)), 100),
                                    number_of_threads=5,
                                    time_cutoff=1000000,
                                    step_cutoff=1000000,
                                    folder=rnsd.param_folder)

print("run simulation with ../RNMC " + network_folder + " " + param_folder)
