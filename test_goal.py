import sys
sys.path.append('./mrnet/src')
from mrnet.network.reaction_generation import *
from monty.serialization import loadfn
import pickle


# you are probably gonna want to do something else while waiting for this
# it takes ~ 30 minutes when running on sams molecule list

molecule_list_json = sys.argv[1]
molecule_entries = loadfn(molecule_list_json)
reaction_generator = ReactionGenerator(molecule_entries)

# don't generate the elementry reactions
reaction_generator.current_chunk = []


with open('./goal','rb') as f:
    goal = pickle.load(f)

missing_reactions = set(goal)


# suggesting that interpreter garbage collect the large frozen set
# trying to keep memory footprint low
goal = None
extras = []
mass_not_conserved = []

def mass_balancer(reaction):
    reactant_atoms = {}
    product_atoms = {}
    for reactant in reaction.reactants:
        for atom in reactant.species:
            if atom in reactant_atoms:
                reactant_atoms[atom] += 1
            else:
                reactant_atoms[atom] = 1

    for product in reaction.products:
        for atom in product.species:
            if atom in product_atoms:
                product_atoms[atom] += 1
            else:
                product_atoms[atom] = 1

    return reactant_atoms, product_atoms



for reaction in reaction_generator:
    reaction_sig = ((frozenset(reaction.reactant_indices),
                     frozenset(reaction.product_indices)))

    if reaction_sig in missing_reactions:
        missing_reactions.remove(reaction_sig)
    else:
        extras.append(reaction_sig)

    reactant_atoms, product_atoms = mass_balancer(reaction)
    if reactant_atoms != product_atoms:
        mass_not_conserved.append(reaction)



if (len(missing_reactions) == 0 and len(extras) == 0 and len(mass_not_conserved) == 0):
    print("all good!")
else:
    print(len(missing_reactions),
          " reactions are missing.")
    print(len(extras),
          " extra reactions.")
    print(len(mass_not_conserved),
          " reactions where atom numbers are not conserved")
    print("this is bad....")
