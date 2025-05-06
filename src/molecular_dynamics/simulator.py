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
from .binding_site import BindingSite


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
            "Usage: python simulator.py <simulate|minimize> <pdb_path> <simulation_params_file> [neighbourhood_radius] [neighbourhood_depth] [chains_file]"
        )
        exit(1)

    command = sys.argv[1]
    if command not in ["simulate", "minimize"]:
        print("Invalid command. Use 'simulate' or 'minimize'.")
        return

    pdb_path = sys.argv[2]
    cmd.load(pdb_path)
    molecule_name = Path(pdb_path).stem

    simulation_params_file = sys.argv[3]
    params = SimulationParameters()
    with open(simulation_params_file, "r") as f:
        yaml_str = f.read()
        params.parse_yaml(yaml_str)

    if command == "minimize":
        simulation = AllAtomSimulationHandler(Path(pdb_path).parent, params)
        fixed_molecule = f"{molecule_name}_fixed"
        simulation.fix_pdb(molecule_name, fixed_molecule)
        simulation.minimize(fixed_molecule, f"{molecule_name}_minimized")
        return

    try:
        neighbourhood_radius = int(sys.argv[4])
        neighbourhood_depth = int(sys.argv[5])

        try:
            chains_file = sys.argv[6]
        except IndexError:
            print("Missing heavy and light chain information.")
            return

        with open(chains_file, "r") as f:
            data = json.load(f)
            heavy_chains = data["heavy_chains"]
            light_chains = data["light_chains"]

        tmp_dir = setup_tmp_dir(
            pdb_path, params, f"r{neighbourhood_radius}d{neighbourhood_depth}"
        )
        simulation = AllAtomSimulationHandler(tmp_dir, params)
        binding_site = BindingSite(molecule_name, heavy_chains, light_chains)
        binding_site.select_paratope()
        binding_site.select_epitope()
        binding_site.update_neighbourhoods(neighbourhood_radius, neighbourhood_depth)
        simulation.simulate_partial(molecule_name, binding_site)

    except IndexError:
        tmp_dir = setup_tmp_dir(pdb_path, params, "full")
        simulation = AllAtomSimulationHandler(tmp_dir, params)
        fixed_molecule = f"{molecule_name}_fixed"
        simulation.fix_pdb(molecule_name, fixed_molecule)
        simulation.simulate(fixed_molecule)

    print(f"Done! Simulation files saved at {tmp_dir}")


if __name__ == "__main__":
    main()
