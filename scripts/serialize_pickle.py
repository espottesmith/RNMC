import os
import sys
import pickle
import random
from pathlib import Path

from pymatgen.core.structure import Molecule
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender
from pymatgen.analysis.graphs import MoleculeGraph

from mrnet.network.reaction_network import *
from rnmc.rnmc import SerializedReactionNetwork, serialize_simulation_parameters


if len(sys.argv) != 4:
    print("usage: python serialize.py pickle network_folder param_folder")
    exit()

network_pickle = sys.argv[1]
network_folder = sys.argv[2]
param_folder = sys.argv[3]

reaction_network = None
with open(network_pickle, "rb") as pick:
    reaction_network = pickle.load(pick)

initial_state = [
    ('../data/molecules/Li.xyz', 1, 600),
    ('../data/molecules/EC.xyz', 0, 7850),
    ('../data/molecules/H2O.xyz', 0, 10),
    ('../data/molecules/H.xyz', 1, 10),
    ('../data/molecules/OH.xyz', -1, 10),
    ('../data/molecules/CO2.xyz', 0, 50)
    ]

initial_state_data = list()
for (mol, charge, number) in initial_state:
    mg = metal_edge_extender(MoleculeGraph.with_local_env_strategy(Molecule.from_file(mol), OpenBabelNN()))
    mg.molecule.set_charge_and_spin(charge)
    initial_state_data.append((mg, number))

rnsd = SerializedReactionNetwork(reaction_network,
                                 initial_state_data,
                                 Path(network_folder),
                                 Path(param_folder),
                                 logging=True)
rnsd.serialize()

if not os.path.exists(rnsd.param_folder):
    serialize_simulation_parameters(seeds=random.sample(list(range(1, 1000000000)), 25),
                                    number_of_threads=6,
                                    time_cutoff=1000000,
                                    step_cutoff=5000000,
                                    folder=rnsd.param_folder)

print("run simulation with ../RNMC " + network_folder + " " + param_folder)
