#!/usr/bin/env python

import sys
import pickle
import os
import argparse
import textwrap
from pathlib import Path

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender

from rnmc.rnmc import SimulationAnalyser


def main():
    parser = argparse.ArgumentParser(
        description="Identify pathways to a particular product of interest",
        epilog=textwrap.dedent("""\
        Example of use:
        ---------------
        search_pathways.py <network_folder> <molecule_file> -c <charge>
        """)
    )
    parser.add_argument("network_folder", help="Folder holding simulation data")
    parser.add_argument("molecule_file", help="Structure file (*.xyz, *.mol, etc.) for target molecule")
    parser.add_argument("-c", "--charge", help="Charge of target molecule", default=0, required=False)

    args = parser.parse_args()

    network_folder = Path(args.network_folder).resolve()

    with open((network_folder / "rnsd.pickle").as_posix(), "rb") as f:
        rnsd = pickle.load(f)

    analyzer = SimulationAnalyser(rnsd, Path(args.network_folder).resolve())

    mg = metal_edge_extender(
        MoleculeGraph.with_local_env_strategy(
            Molecule.from_file(args.molecule_file),
            OpenBabelNN()
        )
    )
    mg.molecule.set_charge_and_spin(int(args.charge))

    index = rnsd.find_index_from_mol_graph(mg)

    if index is None:
        print("TARGET COULD NOT BE FOUND")

    analyzer.generate_pathway_report(index)
