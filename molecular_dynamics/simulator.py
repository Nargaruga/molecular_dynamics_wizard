#!/usr/bin/env python3

import os
import shutil
import time
import sys
import json

from pymol import cmd

from .simulation_params import SimulationParameters
from .aa_simulation_handler import AllAtomSimulationHandler


def simulate_partial_molecule(
    tmp_dir, molecule_name, params, residues_to_simulate, residues_to_lock
):
    simulation = AllAtomSimulationHandler(tmp_dir, params)
    cmd.load(os.path.join("inputs", f"{molecule_name}.pdb"))

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


def setup_tmp_dir(molecule, params, extra=""):
    tmp_dir_path = os.path.join(
        "tmp",
        f"{molecule}_{time.strftime('%Y-%m-%d_%H-%M-%S')}_s{params.sim_steps}_{extra}",
    )
    os.makedirs(tmp_dir_path)
    shutil.copy(f"inputs/{molecule}.pdb", os.path.join(tmp_dir_path, f"{molecule}.pdb"))

    return tmp_dir_path


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python simulator.py <molecule_name> <simulation_params_file> [constraints_file]"
        )
        exit(1)

    molecule = sys.argv[1]

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

        tmp_dir = setup_tmp_dir(molecule, params, f"d{depth}")
        simulate_partial_molecule(
            tmp_dir,
            molecule,
            params,
            set(paratope + paratope_neighbourhood + epitope + epitope_neighbourhood),
            set(locked_paratope_neighbourhood + locked_epitope_neighbourhood),
        )

    except IndexError:
        tmp_dir = setup_tmp_dir(molecule, params, "full")
        simulate_full_molecule(tmp_dir, molecule, params)

    print(f"Done! Simulation files saved at {tmp_dir}")


if __name__ == "__main__":
    main()
