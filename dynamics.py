import os
from enum import Enum, IntEnum, auto
import threading
import time

import json

from pymol.wizard import Wizard
from pymol import cmd

from .molecular_dynamics.aa_simulation_handler import AllAtomSimulationHandler
from .molecular_dynamics.simulation_params import SimulationParameters

# TODO
WIZARD_ROOT = (
    "/home/leo/anaconda3/envs/md_test/lib/python3.9/site-packages/pymol/wizard/"
)


def load_configuration():
    params = SimulationParameters()

    install_data_path = os.path.join("dynamics_extra", "installation_data.json")
    with open(install_data_path) as f:
        data = json.load(f)
        try:
            config_path = data["config_path"]
            params.parse_file(os.path.join(config_path, "simulation_params.yaml"))
        except KeyError:
            print(
                f"WARNING: Failed to read configuration file path from {install_data_path}, using default config."
            )

    return params


class WizardState(IntEnum):
    """The possible states of the wizard."""

    INITIALIZING = auto()
    READY = auto()
    MOLECULE_SELECTED = auto()
    RUNNING_SIMULATION = auto()
    SIMULATION_COMPLETE = auto()


class SimulationType(Enum):
    ALL_ATOM = auto()
    COARSE = auto()


class Residue:
    """A residue in a molecule."""

    def __init__(self, resi, chain):
        self.resi = resi
        self.chain = chain


class Dynamics(Wizard):
    """A wizard for performing molecular dynamics simulations."""

    def __init__(self, _self=cmd):
        Wizard.__init__(self, _self)
        os.chdir(WIZARD_ROOT)

        cmd.set("retain_order", 1)
        cmd.set("pdb_retain_ids", 1)

        self.status = WizardState.INITIALIZING
        self.molecule = None
        self.sim_depth = 2
        self.sim_type = SimulationType.ALL_ATOM
        self.populate_molecule_choices()
        self.populate_sim_depth_choices()
        self.populate_sim_type_choices()
        self.status = WizardState.READY

    def get_prompt(self):
        """Return the prompt for the current state of the wizard."""

        self.prompt = None
        if self.status == WizardState.INITIALIZING:
            self.prompt = ["Initializing, please wait..."]
        elif self.status == WizardState.READY:
            self.prompt = ["Select a molecule."]
        elif self.status == WizardState.MOLECULE_SELECTED:
            self.prompt = [f"Run to perform a simulation for {self.molecule}."]
        elif self.status == WizardState.RUNNING_SIMULATION:
            self.prompt = ["Running simulation, please wait..."]
        elif self.status == WizardState.SIMULATION_COMPLETE:
            self.prompt = ["Simulation complete."]

        return self.prompt

    def get_panel(self):
        """Return the menu panel for the wizard."""

        if self.molecule is None:
            molecule_label = "Choose molecule"
        else:
            molecule_label = self.molecule

        depth_label = f"Neighbourhood: {self.sim_depth}"
        sim_type_label = f"Simulation type: {self.sim_type.name}"

        options = [[1, "Molecular Dynamics", ""], [3, sim_type_label, "sim_type"]]

        if self.status >= WizardState.READY:
            options.extend(
                [[3, molecule_label, "molecule"], [3, depth_label, "sim_depth"]]
            )

        if self.status >= WizardState.MOLECULE_SELECTED:
            options.append(
                [2, "Run Simulation", "cmd.get_wizard().run()"],
            )

        options.append([2, "Dismiss", "cmd.set_wizard()"])

        return options

    def populate_molecule_choices(self):
        """Populate the menu with the available molecules in the session."""

        molecules = cmd.get_names("objects")
        self.menu["molecule"] = [[2, "Molecule", ""]]
        for m in molecules:
            self.menu["molecule"].append(
                [
                    1,
                    m,
                    'cmd.get_wizard().set_molecule("' + m + '")',
                ]
            )

    def populate_sim_depth_choices(self):
        """Populate the menu with the possible values for the depth of the paratope neighbourhood to simulate."""

        self.menu["sim_depth"] = [[2, "Neighbourhood Residues", ""]]
        depths = [0, 2, 4, 8, 10]
        for d in depths:
            self.menu["sim_depth"].append(
                [
                    1,
                    str(d),
                    "cmd.get_wizard().set_sim_depth(" + str(d) + ")",
                ]
            )

    def populate_sim_type_choices(self):
        """Populate the menu with the possible values for the type of simulation to perform."""

        self.menu["sim_type"] = [[2, "Simulation Type", ""]]
        for sim_type in SimulationType:
            self.menu["sim_type"].append(
                [
                    1,
                    sim_type.name,
                    'cmd.get_wizard().set_sim_type("' + sim_type.name + '")',
                ]
            )

    def set_molecule(self, molecule):
        """Set the molecule to simulate."""

        self.molecule = molecule
        self.status = WizardState.MOLECULE_SELECTED
        cmd.refresh_wizard()

    def set_sim_depth(self, depth):
        """Set the depth of the paratope neighbourhood to simulate."""

        self.sim_depth = depth
        cmd.refresh_wizard()

    def set_sim_type(self, sim_type_str):
        """Set the type of simulation to perform."""

        if sim_type_str == "ALL_ATOM":
            self.sim_type = SimulationType.ALL_ATOM
        elif sim_type_str == "COARSE":
            self.sim_type = SimulationType.COARSE
        cmd.refresh_wizard()

    def run(self):
        """Run the wizard to perform a molecular dynamics simulation of the selected molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        self.status = WizardState.RUNNING_SIMULATION
        cmd.refresh_wizard()

        # Run the simulation on a separate thread to keep the interface responsive
        worker_thread = threading.Thread(
            target=self.run_simulation, args=[self.sim_depth]
        )
        worker_thread.start()

    def run_simulation(self, depth=2, reps=1, cleanup=False):
        """Run the simulation."""

        if self.molecule is None:
            print("Please select a molecule.")
            self.status = WizardState.READY
            cmd.refresh_wizard()
            return

        sim_params = load_configuration()

        tmp_dir = os.path.join(
            "tmp",
            f"{self.molecule}_{time.strftime('%Y-%m-%d_%H-%M-%S')}_s{sim_params.sim_steps}_d{depth}",
        )
        os.makedirs(tmp_dir)

        simulation = AllAtomSimulationHandler(tmp_dir, sim_params)

        # Get all residues in the paratope
        try:
            self.identify_paratope(self.molecule, f"{self.molecule}_paratope")
        except Exception as _:
            cmd.set_wizard()
            print("Could not identify the paratope. Is the molecule an antibody?")
            self.status = WizardState.READY
            cmd.refresh_wizard()
            return

        self.select_epitope(
            self.molecule, f"{self.molecule}_paratope", f"{self.molecule}_epitope"
        )

        paratope_residues = set()
        cmd.iterate(
            f"{self.molecule}_paratope",
            "paratope_residues.add((resi, chain))",
            space=locals(),
        )

        epitope_residues = set()
        cmd.iterate(
            f"{self.molecule}_epitope",
            "epitope_residues.add((resi, chain))",
            space=locals(),
        )

        # Get the neighbourhood of the paratope
        paratope_neigh_residues = self.get_neigbourhood(
            f"{self.molecule}_paratope_neigh",
            f"{self.molecule}_paratope",
            "chain H or chain L",
            depth,
        )

        # Obtain a further extended neighbourhood of the paratope that will have its atoms locked
        locked_paratope_neigh_residues = self.get_neigbourhood(
            f"{self.molecule}_locked_paratope_neigh",
            f"{self.molecule}_paratope_neigh",
            "chain H or chain L",
        )

        # Get the neighbourhood of the epitope
        epitope_neigh_residues = self.get_neigbourhood(
            f"{self.molecule}_epitope_neigh",
            f"{self.molecule}_epitope",
            "not chain H and not chain L",
            depth,
        )

        # Obtain a further extended neighbourhood of the epitope that will have its atoms locked
        locked_epitope_neigh_residues = self.get_neigbourhood(
            f"{self.molecule}_locked_epitope_neigh",
            f"{self.molecule}_epitope_neigh",
            "not chain H and not chain L",
        )

        # Create a new object out of the non-simulated atoms
        non_sim_molecule = f"{self.molecule}_non_sim"
        non_sim_molecule_sel = f"{non_sim_molecule}_sel"
        # This object includes the locked neighbourhood, which will help us
        # align the final trajaectory with the original structure
        cmd.select(
            name=non_sim_molecule_sel,
            selection=f"byres {self.molecule} and not {self.molecule}_paratope_neigh and not {self.molecule}_epitope_neigh",
        )
        cmd.create(non_sim_molecule, non_sim_molecule_sel)
        cmd.delete(non_sim_molecule_sel)
        cmd.align(non_sim_molecule, self.molecule)
        cmd.disable(non_sim_molecule)

        # Save residue ids to json
        with open(os.path.join(tmp_dir, "simulated_residues.json"), "w") as f:
            json.dump(
                {
                    "paratope": list(paratope_residues),
                    "paratope_neighbourhood": list(paratope_neigh_residues),
                    "locked_paratope_neighbourhood": list(
                        locked_paratope_neigh_residues
                    ),
                    "epitope": list(epitope_residues),
                    "epitope_neighbourhood": list(epitope_neigh_residues),
                    "locked_epitope_neighbourhood": list(locked_epitope_neigh_residues),
                    "depth": depth,
                },
                f,
            )

        if sim_params.remove_non_simulated:
            print("Removing non-simulated atoms...")
            # This operation causes the renumbering of all atoms,
            # but the residue ids are preserved
            simulation.slice_object(
                self.molecule,
                locked_paratope_neigh_residues.union(locked_epitope_neigh_residues),
            )

            to_fix = f"{self.molecule}_sliced"
        else:
            to_fix = self.molecule

        cmd.save(os.path.join(tmp_dir, f"{to_fix}.pdb"), to_fix)

        final_molecule = f"{self.molecule}_fixed"
        simulation.preprocess_input(to_fix, f"{final_molecule}.pdb")
        cmd.disable(self.molecule)
        cmd.disable(to_fix)
        cmd.load(os.path.join(tmp_dir, f"{final_molecule}.pdb"))
        cmd.show_as("licorice", final_molecule)

        # Remove overlap
        locked_paratope_neigh_residues = (
            locked_paratope_neigh_residues - paratope_neigh_residues
        )
        paratope_neigh_residues = paratope_neigh_residues - paratope_residues

        locked_epitope_neigh_residues = (
            locked_epitope_neigh_residues - epitope_neigh_residues
        )
        epitope_neigh_residues = epitope_neigh_residues - epitope_residues

        # Run the simulation
        simulation.simulate(
            final_molecule,
            simulation.residues_to_atoms(
                final_molecule,
                paratope_residues.union(paratope_neigh_residues)
                .union(epitope_residues)
                .union(epitope_neigh_residues),
            ),
        )
        simulation.postprocess_output()

        output_name = final_molecule
        cmd.load_traj(os.path.join(tmp_dir, "trajectory.dcd"), output_name)

        # Colour the residues
        cmd.color("grey", f"{output_name}")

        for resi, chain in locked_paratope_neigh_residues:
            cmd.color("red", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in paratope_neigh_residues:
            cmd.color("blue", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in paratope_residues:
            cmd.color("green", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in locked_epitope_neigh_residues:
            cmd.color("magenta", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in epitope_neigh_residues:
            cmd.color("cyan", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in epitope_residues:
            cmd.color("yellow", f"{output_name} and resi {resi} and chain {chain}")

        cmd.align(output_name, self.molecule)
        cmd.zoom(output_name)

        if cleanup:
            print("Cleaning up...")

            # Remove molecules
            cmd.delete(self.molecule)
            cmd.delete(self.molecule)
            cmd.delete(non_sim_molecule)

            # Remove selections
            cmd.delete(f"{self.molecule}_paratope")
            cmd.delete(f"{self.molecule}_paratope_neigh")
            cmd.delete(f"{self.molecule}_locked_paratope_neigh")

        self.status = WizardState.SIMULATION_COMPLETE
        cmd.refresh_wizard()

        print(f"Done! Simulation files saved at {tmp_dir}")

    def get_neigbourhood(self, selection_name, target_name, chains, depth=1):
        residues = set()

        if depth == 0:
            cmd.select(name=selection_name, selection=target_name)
        else:
            current_selection = target_name

            for _ in range(depth):
                cmd.select(
                    name=selection_name,
                    selection=f"byres ((neighbor {current_selection}) or {current_selection}) and ({chains})",
                    merge=1,
                )

                current_selection = selection_name

        cmd.iterate(
            selection_name,
            "residues.add((resi, chain))",
            space=locals(),
        )

        return residues

    def identify_paratope(self, molecule, selection_name):
        """Identify the paratope on the selected antibody through the appropriate wizard."""

        cmd.wizard("paratope")
        cmd.get_wizard().set_molecule(molecule)
        cmd.get_wizard().set_selection_name(selection_name)
        cmd.get_wizard().run()
        cmd.get_wizard().toggle_label_pos()
        cmd.set_wizard()

    def select_epitope(self, molecule, paratope_sel, selection_name):
        """Get the residues of the epitope based on the paratope selection."""

        epitope_sel = selection_name
        cmd.select(
            name=epitope_sel,
            selection=f"{molecule} and not chain H and not chain L near_to 6.0 of {paratope_sel}",
        )
