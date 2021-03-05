import sys
sys.path.append('./mrnet/src')
from mrnet.network.reaction_network import *
from monty.serialization import loadfn
import pickle

molecule_list_json = sys.argv[1]
molecule_entries = loadfn(molecule_list_json)

rn = ReactionNetwork.from_input_entries(molecule_entries)
rn.build()
pre_goal = ReactionNetwork.identify_concerted_rxns_via_intermediates(
    rn,[e.parameters["ind"] for e in rn.entries_list])

def filter_nones(list):
    result = []
    for x in list:
        if x != None:
            result.append(x)

    return result

goal_list = []

for (reactants, products, _) in pre_goal:
    sig = (frozenset(filter_nones(reactants)), frozenset(filter_nones(products)))
    goal_list.append(sig)


goal = frozenset(goal_list)

with open('./goal','wb') as f:
    pickle.dump(goal,f)
