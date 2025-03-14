#!/usr/bin/env python3

import os
from pathlib import Path
import shutil
import time
import sys
import json

from pymol import cmd

from .simulation_params import SimulationParameters
from .aa_simulation_handler import AllAtomSimulationHandler


def simulate_partial_molecule(
    tmp_dir, pdb_path, params, residues_to_simulate, residues_to_lock
):
    simulation = AllAtomSimulationHandler(tmp_dir, params)
    cmd.load(pdb_path)

    molecule_name = Path(pdb_path).stem
    to_fix = molecule_name
    if params.remove_non_simulated:
        residues_to_keep = residues_to_simulate.union(residues_to_lock)
        simulation.slice_object(molecule_name, residues_to_keep)
        to_fix = f"{molecule_name}_sliced"

    cmd.save(os.path.join(tmp_dir, f"{to_fix}.pdb"), to_fix)
    final_molecule = f"{molecule_name}_fixed"
    simulation.preprocess_input(to_fix, f"{final_molecule}.pdb")
    cmd.load(os.path.join(tmp_dir, f"{final_molecule}.pdb"))

    simulation.simulate(
        final_molecule,
        simulation.residues_to_atoms(final_molecule, residues_to_simulate),
    )
    simulation.postprocess_output()


def simulate_full_molecule(tmp_dir, molecule_name, params):
    simulation = AllAtomSimulationHandler(tmp_dir, params)

    simulation.preprocess_input(molecule_name, f"{molecule_name}_fixed.pdb")
    simulation.simulate(f"{molecule_name}_fixed")
    simulation.postprocess_output()


def lists_to_tuples(lists):
    return [tuple(lst) for lst in lists]


def setup_tmp_dir(pdb_path, params, extra=""):
    tmp_dir_path = os.path.join(
        "results",
        f"{Path(pdb_path).stem}_{time.strftime('%Y-%m-%d_%H-%M-%S')}_s{params.sim_steps}_{extra}",
    )
    os.makedirs(tmp_dir_path)
    shutil.copy(pdb_path, os.path.join(tmp_dir_path, Path(pdb_path).name))

    return tmp_dir_path


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python simulator.py <pdb_path> <simulation_params_file> [constraints_file]"
        )
        exit(1)

    pdb_path = sys.argv[1]

    simulation_params_file = sys.argv[2]
    params = SimulationParameters()
    params.parse_file(simulation_params_file)

    try:
        constraints_file = sys.argv[3]
        with open(constraints_file, "r") as file:
            data = json.load(file)
            paratope = lists_to_tuples(data["paratope"])
            paratope_neighbourhood = lists_to_tuples(data["paratope_neighbourhood"])
            locked_paratope_neighbourhood = lists_to_tuples(
                data["locked_paratope_neighbourhood"]
            )
            epitope = lists_to_tuples(data["epitope"])
            epitope_neighbourhood = lists_to_tuples(data["epitope_neighbourhood"])
            locked_epitope_neighbourhood = lists_to_tuples(
                data["locked_epitope_neighbourhood"]
            )
            depth = data["depth"]

        tmp_dir = setup_tmp_dir(params, f"d{depth}")
        simulate_partial_molecule(
            tmp_dir,
            pdb_path,
            params,
            set(paratope + paratope_neighbourhood + epitope + epitope_neighbourhood),
            set(locked_paratope_neighbourhood + locked_epitope_neighbourhood),
        )

    except IndexError:
        tmp_dir = setup_tmp_dir(pdb_path, params, "full")
        simulate_full_molecule(tmp_dir, Path(pdb_path).stem, params)

    print(f"Done! Simulation files saved at {tmp_dir}")


if __name__ == "__main__":
    main()
