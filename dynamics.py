import os
from enum import Enum, IntEnum, auto
import threading
import time
import pathlib

from pymol.wizard import Wizard
from pymol import cmd

from .molecular_dynamics.aa_simulation_handler import AllAtomSimulationHandler
from .molecular_dynamics.simulation_params import SimulationParameters


def load_configuration(installed_wizard_path):
    params = SimulationParameters()

    config_path = os.path.join(
        installed_wizard_path, "dynamics_extra", "simulation_params.yaml"
    )
    params.parse_file(config_path)

    return params


def get_installed_wizard_path():
    return pathlib.Path(__file__).parent.resolve()


class WizardState(IntEnum):
    """The possible states of the wizard."""

    INITIALIZING = auto()
    READY = auto()
    MOLECULE_SELECTED = auto()
    CHAINS_SELECTED = auto()
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
        os.chdir(get_installed_wizard_path())

        cmd.set("retain_order", 1)
        cmd.set("pdb_retain_ids", 1)

        self.status = WizardState.INITIALIZING
        self.molecule = None
        self.heavy_chains = []
        self.light_chains = []
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
            self.prompt = [f"Select the heavy and light chains for {self.molecule}."]
        elif self.status == WizardState.CHAINS_SELECTED:
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

        heavy_chain_label = "Heavy Chains: "
        if self.heavy_chains:
            heavy_chain_label += ", ".join(self.heavy_chains)
        else:
            heavy_chain_label += "None"

        light_chain_label = "Light Chains: "
        if self.light_chains:
            light_chain_label += ", ".join(self.light_chains)
        else:
            light_chain_label += "None"

        if self.status >= WizardState.MOLECULE_SELECTED:
            options.extend(
                [
                    [3, heavy_chain_label, "heavy_chain"],
                    [3, light_chain_label, "light_chain"],
                ]
            )

        if self.status >= WizardState.CHAINS_SELECTED:
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

    def populate_chain_choices(self):
        """Populate the menu with the available chains in the selected molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        chains = cmd.get_chains(self.molecule)
        self.menu["heavy_chain"] = [[2, "Heavy Chain", ""]]
        for c in chains:
            self.menu["heavy_chain"].append(
                [
                    1,
                    c,
                    'cmd.get_wizard().set_heavy_chain("' + c + '")',
                ]
            )

        self.menu["light_chain"] = [[2, "Light Chain", ""]]
        for c in chains:
            self.menu["light_chain"].append(
                [
                    1,
                    c,
                    'cmd.get_wizard().set_light_chain("' + c + '")',
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
        self.populate_chain_choices()
        cmd.refresh_wizard()

    def set_heavy_chain(self, chain):
        """Set the heavy chain to be used for the heatmap."""

        if chain in self.heavy_chains:
            self.heavy_chains.remove(chain)
        else:
            self.heavy_chains.append(chain)

        if self.heavy_chains and self.light_chains:
            self.status = WizardState.CHAINS_SELECTED

        cmd.refresh_wizard()

    def set_light_chain(self, chain):
        """Set the light chain to be used for the heatmap."""

        if chain in self.light_chains:
            self.light_chains.remove(chain)
        else:
            self.light_chains.append(chain)

        if self.heavy_chains and self.light_chains:
            self.status = WizardState.CHAINS_SELECTED

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

    def run_simulation(self, depth=2):
        """Run the simulation."""

        if self.molecule is None:
            print("Please select a molecule.")
            self.status = WizardState.READY
            cmd.refresh_wizard()
            return

        sim_params = load_configuration(get_installed_wizard_path())
        tmp_dir = os.path.join(
            "simulations",
            f"{self.molecule}_{time.strftime('%Y-%m-%d_%H-%M-%S')}_s{sim_params.sim_steps}_d{depth}",
        )
        os.makedirs(tmp_dir)
        cmd.save(os.path.join(tmp_dir, f"{self.molecule}.pdb"), self.molecule)
        try:
            simulation = AllAtomSimulationHandler(tmp_dir, sim_params)
            simulation.simulate(
                self.molecule, depth, self.heavy_chains, self.light_chains
            )
        except Exception as e:
            cmd.set_wizard()
            print(f"Error while running simulation: {e}.")
            self.status = WizardState.READY
            cmd.refresh_wizard()
            return

        simulation_data = simulation.simulation_data
        constrained_residues = set()
        for atom in simulation_data.constrained_atoms:
            cmd.select(
                "constrained_atoms",
                f"byres {simulation.simulation_data.final_molecule} and index {atom}",
            )
        cmd.iterate(
            "constrained_atoms",
            "constrained_residues.add((resi, chain))",
            space=locals(),
        )

        output_name = simulation.simulation_data.final_molecule
        cmd.load_traj(os.path.join(tmp_dir, "trajectory.dcd"), output_name)

        # Colour the residues
        cmd.color("grey", f"{output_name}")

        for resi, chain in simulation_data.locked_paratope_neigh_residues:
            cmd.color("red", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in simulation_data.paratope_neigh_residues:
            cmd.color("blue", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in simulation_data.paratope_residues:
            cmd.color("green", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in simulation_data.locked_epitope_neigh_residues:
            cmd.color("magenta", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in simulation_data.epitope_neigh_residues:
            cmd.color("cyan", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in simulation_data.epitope_residues:
            cmd.color("yellow", f"{output_name} and resi {resi} and chain {chain}")

        for resi, chain in constrained_residues:
            cmd.color("black", f"{output_name} and resi {resi} and chain {chain}")

        cmd.align(output_name, self.molecule)
        cmd.zoom(output_name)

        self.status = WizardState.SIMULATION_COMPLETE
        cmd.refresh_wizard()

        print(f"Done! Simulation files saved at {tmp_dir}")
