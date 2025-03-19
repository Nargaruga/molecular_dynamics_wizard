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
            "Usage: python simulator.py <pdb_path> <simulation_params_file> [neighbourhood_depth] [chains_file]"
        )
        exit(1)

    pdb_path = sys.argv[1]
    cmd.load(pdb_path)
    molecule_name = Path(pdb_path).stem

    simulation_params_file = sys.argv[2]
    params = SimulationParameters()
    params.parse_file(simulation_params_file)

    try:
        neighbourhood_depth = sys.argv[3]

        try:
            chains_file = sys.argv[4]
        except IndexError:
            print("Missing heavy and light chain information.")
            return

        with open(chains_file, "r") as f:
            data = json.load(f)
            heavy_chains = data["heavy_chains"]
            light_chains = data["light_chains"]

        tmp_dir = setup_tmp_dir(params, f"d{neighbourhood_depth}")
        simulation = AllAtomSimulationHandler(tmp_dir, params)
        simulation.simulate_partial(
            molecule_name, neighbourhood_depth, heavy_chains, light_chains
        )

    except IndexError:
        tmp_dir = setup_tmp_dir(pdb_path, params, "full")
        simulation = AllAtomSimulationHandler(tmp_dir, params)
        simulation.simulate(molecule_name)

    print(f"Done! Simulation files saved at {tmp_dir}")


if __name__ == "__main__":
    main()
