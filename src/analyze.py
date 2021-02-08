from rnmc.rnmc import *
import sys
import pickle


if len(sys.argv) != 2:
    print("usage: python analyze.py network_folder")
    exit()

network_folder = sys.argv[1]

with open(network_folder + "/rnsd.pickle",'rb') as f:
    rnsd = pickle.load(f)

# TODO: don't want network folder as an attribute of rnsd
rnsd.network_folder = network_folder
sa = SimulationAnalyser(rnsd, network_folder)


LEDC_index = rnsd.find_index_from_mol_graph('./mrnet/test_files/reaction_network_files/LEDC.xyz',0)
Li2CO3_index = rnsd.find_index_from_mol_graph('./li2co3_0.xyz', 0)
sa.generate_pathway_report(LEDC_index)
