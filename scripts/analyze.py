import sys
import pickle
import os

from rnmc.rnmc import SimulationAnalyser


if len(sys.argv) != 2:
    print("usage: python analyze.py network_folder target charge")
    exit()

network_folder = sys.argv[1]

with open(os.path.join(network_folder, "rnsd.pickle"), "rb") as f:
    rnsd = pickle.load(f)

# TODO: don't want network folder as an attribute of rnsd
rnsd.network_folder = network_folder
sa = SimulationAnalyser(rnsd, network_folder)

target = sys.argv[2]
charge = sys.argv[3]

index = rnsd.find_index_from_mol_graph(target, int(charge))
sa.generate_pathway_report(index)
