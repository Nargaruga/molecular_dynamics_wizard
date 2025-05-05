import os
from enum import Enum, IntEnum, auto
import threading
import time
import tempfile
import pkgutil
import threading

from pymol.wizard import Wizard
from pymol.wizarding import WizardError
from pymol import cmd, CmdException

from molecular_dynamics.aa_simulation_handler import AllAtomSimulationHandler
from molecular_dynamics.simulation_params import SimulationParameters
from molecular_dynamics.binding_site import (
    BindingSite,
    BindingSiteError,
    get_residues,
)


def load_configuration() -> SimulationParameters:
    """Load simulation parameters from file and return them."""

    params = SimulationParameters()

    raw_yaml = pkgutil.get_data(
        "molecular_dynamics", os.path.join("config", "simulation_params.yaml")
    )
    if raw_yaml is None:
        raise WizardError("Could not load simulation parameters.")

    yaml_str = raw_yaml.decode("utf-8")
    params.parse_yaml(yaml_str)

    return params


class WizardInputState(IntEnum):
    READY = auto()
    MOLECULE_SELECTED = auto()
    CHAINS_SELECTED = auto()
    PARATOPE_SELECTED = auto()
    SIMULATION_READY = auto()


class WizardTask(IntEnum):
    IDENTIFYING_BINDING_SITE = auto()
    MINIMIZING_ENERGY = auto()
    RUNNING_SIMULATION = auto()


class ParatopeDetectionType(Enum):
    EXISTING = auto()
    NEW = auto()


class SimulationType(Enum):
    FULL = "Full"
    PARTIAL = "Partial"


class Residue:
    """A residue in a molecule."""

    def __init__(self, resi, chain):
        self.resi = resi
        self.chain = chain


class Dynamics(Wizard):
    """A wizard for performing molecular dynamics simulations."""

    def __init__(self, _self=cmd):
        Wizard.__init__(self, _self)

        cmd.set("retain_order", 1)
        cmd.set("pdb_retain_ids", 1)

        self.sim_params = load_configuration()
        self.molecule = None
        self.heavy_chains = []
        self.light_chains = []
        self.sim_radius = 2
        self.sim_depth = 2
        self.sim_type = SimulationType.FULL
        self.populate_molecule_choices()
        self.populate_sim_type_choices()
        self.paratope_detection_type = ParatopeDetectionType.NEW
        self.custom_paratope = None
        self.binding_site = None
        self.input_state = WizardInputState.READY
        self.state_lock = threading.Lock()
        self.tasks = []
        self.tasks_lock = threading.Lock()
        self.extra_msg = ""

    def get_prompt(self):  # type: ignore
        """Return the prompt for the current state of the wizard."""

        prompt = []
        if self.input_state == WizardInputState.READY:
            prompt.append("Select a molecule to continue.")
        elif self.input_state == WizardInputState.MOLECULE_SELECTED:
            prompt.append("Select the heavy and light chains to continue.")
        elif self.input_state == WizardInputState.CHAINS_SELECTED:
            prompt.append("Identify the binding site to continue.")
        elif self.input_state == WizardInputState.SIMULATION_READY:
            prompt.append("The simulation is ready to be started.")

        for task in self.tasks:
            if task == WizardTask.IDENTIFYING_BINDING_SITE:
                prompt.append("Detecting binding site, please wait...")
            elif task == WizardTask.MINIMIZING_ENERGY:
                prompt.append("Minimizing energy, please wait...")
            elif task == WizardTask.RUNNING_SIMULATION:
                prompt.append("Running simulation, please wait...")

        return prompt

    def get_panel(self):  # type: ignore
        # Title
        options = [
            [1, "Molecular Dynamics", ""],
        ]

        # Basic entries
        if self.input_state >= WizardInputState.READY:
            if self.molecule is None:
                molecule_label = "Choose molecule"
            else:
                molecule_label = self.molecule

            sim_type_label = f"Simulation type: {self.sim_type.value}"

            options.extend(
                [
                    [3, molecule_label, "molecule"],
                    [3, sim_type_label, "sim_type"],
                ]
            )

        if self.input_state >= WizardInputState.MOLECULE_SELECTED:
            options.extend(
                [
                    [2, "Minimize Energy", "cmd.get_wizard().minimize_structure()"],
                ]
            )

        # Add entries for partial simulations
        if self.sim_type == SimulationType.PARTIAL:
            if self.input_state >= WizardInputState.MOLECULE_SELECTED:
                paratope_detection_type_label = "Mode: "
                if self.paratope_detection_type == ParatopeDetectionType.NEW:
                    paratope_detection_type_label += "New paratope"
                else:
                    paratope_detection_type_label += "Existing paratope"

                options.append(
                    [
                        2,
                        paratope_detection_type_label,
                        "cmd.get_wizard().toggle_paratope_detection_mode()",
                    ],
                )
                heavy_chains_label = "Heavy Chains: "
                if self.heavy_chains:
                    heavy_chains_label += ", ".join(self.heavy_chains)
                else:
                    heavy_chains_label += "None"

                light_chains_label = "Light Chains: "
                if self.light_chains:
                    light_chains_label += ", ".join(self.light_chains)
                else:
                    light_chains_label += "None"

                options.extend(
                    [
                        [3, heavy_chains_label, "heavy_chain"],
                        [3, light_chains_label, "light_chain"],
                    ]
                )

            detect_binging_site_btn = [
                2,
                "Detect Binding Site",
                "cmd.get_wizard().detect_binding_site()",
            ]

            if self.input_state >= WizardInputState.CHAINS_SELECTED:
                if self.paratope_detection_type == ParatopeDetectionType.EXISTING:
                    paratope_selection_label = "Paratope: "
                    if self.custom_paratope:
                        paratope_selection_label += self.custom_paratope
                    else:
                        paratope_selection_label += "None"
                    options.append(
                        [3, paratope_selection_label, "paratope_selection"],
                    )
                elif self.paratope_detection_type == ParatopeDetectionType.NEW:
                    options.append(detect_binging_site_btn)

            if self.input_state >= WizardInputState.PARATOPE_SELECTED:
                options.append(detect_binging_site_btn)

            if self.input_state >= WizardInputState.SIMULATION_READY:
                radius_label = f"Neighbourhood Radius: {self.sim_radius}"
                depth_label = f"Neighbourhood Depth: {self.sim_depth}"

                options.extend(
                    [
                        [3, radius_label, "sim_radius"],
                        [3, depth_label, "sim_depth"],
                    ]
                )

        if self.input_state >= WizardInputState.SIMULATION_READY:
            options.append(
                [2, "Run Simulation", "cmd.get_wizard().run()"],
            )

        options.append([2, "Dismiss", "cmd.set_wizard()"])

        return options

    def add_task(self, task: WizardTask):
        """Add a task to the list of tasks."""

        if self.check_task(task):
            return

        with self.tasks_lock:
            self.tasks.append(task)
            cmd.refresh_wizard()

    def remove_task(self, task: WizardTask):
        """Remove a task from the list of tasks."""

        if not self.check_task(task):
            return

        with self.tasks_lock:
            self.tasks.remove(task)
            cmd.refresh_wizard()

    def check_task(self, task: WizardTask) -> bool:
        """Check if a task is in the list of tasks."""

        with self.tasks_lock:
            return task in self.tasks

    def update_input_state(self):
        """Update the state of the wizard based on the current inputs."""

        if self.molecule is None:
            self.input_state = WizardInputState.READY
            return

        if self.sim_type == SimulationType.FULL:
            # we only need the molecule for a full simulation
            self.input_state = WizardInputState.SIMULATION_READY
        else:
            # partial simulations require additional inputs
            if self.molecule:
                self.input_state = WizardInputState.MOLECULE_SELECTED

            if self.heavy_chains and self.light_chains:
                self.input_state = WizardInputState.CHAINS_SELECTED

            if (
                self.paratope_detection_type == ParatopeDetectionType.EXISTING
                and self.custom_paratope
            ):
                self.input_state = WizardInputState.PARATOPE_SELECTED

            if self.binding_site:
                self.input_state = WizardInputState.SIMULATION_READY

        cmd.refresh_wizard()

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

    def populate_sim_radius_choices(self):
        """Populate the menu with the possible values for the radius of the paratope neighbourhood to simulate."""

        self.menu["sim_radius"] = [[2, "Neighbourhood Radius", ""]]
        radii = [
            0,
            2,
            4,
            8,
        ]  # Angstrom TODO: generate based on the cutoff distance in the simulation parameters
        for r in radii:
            self.menu["sim_radius"].append(
                [
                    1,
                    str(r),
                    "cmd.get_wizard().set_sim_radius(" + str(r) + ")",
                ]
            )

    def populate_sim_depth_choices(self):
        """Populate the menu with the possible values for the depth of the paratope neighbourhood to simulate."""

        self.menu["sim_depth"] = [[2, "Neighbourhood Residues", ""]]
        depths = [
            0,
            2,
            4,
            8,
        ]
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
                    sim_type.value,
                    'cmd.get_wizard().set_sim_type("' + sim_type.value + '")',
                ]
            )

    def populate_paratope_selection_choices(self):
        self.menu["paratope_selection"] = [[2, "Paratope Selection", ""]]
        for selection in cmd.get_names("selections"):
            self.menu["paratope_selection"].append(
                [
                    1,
                    selection,
                    'cmd.get_wizard().set_paratope_selection("' + selection + '")',
                ]
            )

    def set_paratope_selection(self, selection):
        self.custom_paratope = selection
        self.update_input_state()

    def toggle_paratope_detection_mode(self):
        """Toggle between existing and new paratope detection modes."""

        if self.paratope_detection_type == ParatopeDetectionType.NEW:
            self.paratope_detection_type = ParatopeDetectionType.EXISTING
            self.populate_paratope_selection_choices()
        else:
            self.paratope_detection_type = ParatopeDetectionType.NEW
            self.populate_chain_choices()

        self.update_input_state()

    def attempt_autocomplete(self):
        cmd.wizard("paratope")
        cmd.get_wizard().set_molecule(self.molecule)
        chains = cmd.get_wizard().fetch_cached_chains()
        cmd.set_wizard()

        if chains is None:
            return

        if chains.heavy_chains:
            self.heavy_chains = chains.heavy_chains
        if chains.light_chains:
            self.light_chains = chains.light_chains
        self.update_input_state()

    def set_molecule(self, molecule):
        """Set the molecule to simulate."""

        self.molecule = molecule
        self.binding_site = None
        self.heavy_chains = []
        self.light_chains = []
        self.populate_chain_choices()
        self.attempt_autocomplete()
        self.update_input_state()

    def set_heavy_chain(self, chain):
        """Set the heavy chain to be used for the heatmap."""

        if chain in self.heavy_chains:
            self.heavy_chains.remove(chain)
        else:
            self.heavy_chains.append(chain)

        self.update_input_state()

    def set_light_chain(self, chain):
        """Set the light chain to be used for the heatmap."""

        if chain in self.light_chains:
            self.light_chains.remove(chain)
        else:
            self.light_chains.append(chain)

        self.update_input_state()

    def set_sim_radius(self, radius):
        """Set the radius of the paratope neighbourhood to simulate."""

        self.sim_radius = radius
        self.update_neighbourhoods()
        self.update_coloring(self.molecule)
        cmd.refresh_wizard()

    def set_sim_depth(self, depth):
        """Set the depth of the paratope neighbourhood to simulate."""

        self.sim_depth = depth
        self.update_neighbourhoods()
        self.update_coloring(self.molecule)
        cmd.refresh_wizard()

    def set_sim_type(self, sim_type_str):
        """Set the type of simulation to perform."""

        if sim_type_str == "Full":
            self.sim_type = SimulationType.FULL
        elif sim_type_str == "Partial":
            self.sim_type = SimulationType.PARTIAL
            self.populate_sim_radius_choices()
            self.populate_sim_depth_choices()

        self.update_input_state()

    def update_neighbourhoods(self):
        if self.binding_site is None:
            print("Please perform binding site detection first.")
            return

        self.binding_site.update_neighbourhoods(self.sim_radius, self.sim_depth)

    def update_coloring(self, molecule):
        """Update the coloring of the molecule based on the selected neighbourhood."""

        if self.binding_site is None:
            return

        cmd.color("grey", molecule)

        self.color_residues(
            molecule, get_residues(self.binding_site.paratope_sel), "green"
        )
        self.color_residues(
            molecule, get_residues(self.binding_site.paratope_neigh_sel), "blue"
        )
        self.color_residues(
            molecule, get_residues(self.binding_site.ext_paratope_neigh_sel), "red"
        )

        self.color_residues(
            molecule, get_residues(self.binding_site.epitope_sel), "yellow"
        )
        self.color_residues(
            molecule, get_residues(self.binding_site.epitope_neigh_sel), "cyan"
        )
        self.color_residues(
            molecule, get_residues(self.binding_site.ext_epitope_neigh_sel), "magenta"
        )

    def detect_binding_site(self):
        """Detect the binding site of the selected molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        def aux():
            if self.check_task(WizardTask.IDENTIFYING_BINDING_SITE):
                self.extra_msg = (
                    "Binding site detection is already in progress, please wait..."
                )
                print(self.extra_msg)
                return

            self.add_task(WizardTask.IDENTIFYING_BINDING_SITE)
            try:
                binding_site = BindingSite(
                    self.molecule, self.heavy_chains, self.light_chains
                )
                if (
                    self.paratope_detection_type == ParatopeDetectionType.EXISTING
                    and self.custom_paratope
                ):
                    print(f"Using existing paratope selection: {self.custom_paratope}")
                    binding_site.paratope_sel = self.custom_paratope
                binding_site.select(self.sim_radius, self.sim_depth)
                self.binding_site = binding_site
                self.extra_msg = ""
                self.update_coloring(self.molecule)
                print("Binding site highlighted.")
            except BindingSiteError as e:
                print(f"Error while detecting binding site: {e}.")
                self.binding_site = None
                return
            finally:
                self.remove_task(WizardTask.IDENTIFYING_BINDING_SITE)
                self.extra_msg = ""
                self.update_input_state()

        worked_thread = threading.Thread(
            target=aux,
        )
        worked_thread.start()

    def color_residues(self, molecule, residues, color):
        """Colour the residues in the molecule with the specified colour."""

        for resi, chain in residues:
            cmd.color(color, f"{molecule} and resi {resi} and chain {chain}")

    def minimize_structure(self):
        """Minimize the energy of the selected molecule."""

        def aux():
            if self.molecule is None:
                print("Please select a molecule.")
                return

            if self.check_task(WizardTask.MINIMIZING_ENERGY):
                self.extra_msg = (
                    "Energy minimization is already in progress, please wait..."
                )
                print(self.extra_msg)
                return

            self.add_task(WizardTask.MINIMIZING_ENERGY)
            try:
                with tempfile.TemporaryDirectory() as tmp_dir:
                    cmd.save(
                        os.path.join(tmp_dir, f"{self.molecule}.pdb"), self.molecule
                    )

                    simulation = AllAtomSimulationHandler(tmp_dir, self.sim_params)
                    fixed_molecule = "fixed"
                    simulation.fix_pdb(self.molecule, fixed_molecule)
                    minimized_molecule = "minimized"
                    simulation.minimize(fixed_molecule, minimized_molecule)

                    cmd.load(
                        os.path.join(tmp_dir, f"{minimized_molecule}.pdb"),
                        f"{self.molecule}_min",
                    )

                cmd.disable(self.molecule)

                # TODO locking might be necessary here
                self.populate_molecule_choices()
                self.set_molecule(f"{self.molecule}_min")

                print("Energy minimization complete.")
            # TODO: handle exceptions
            finally:
                self.remove_task(WizardTask.MINIMIZING_ENERGY)
                self.extra_msg = ""
                self.update_input_state()

        worker_thread = threading.Thread(
            target=aux,
        )
        worker_thread.start()

    def run(self):
        """Run the wizard to perform a molecular dynamics simulation of the selected molecule."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        # Run the simulation on a separate thread to keep the interface responsive
        def aux():
            if self.check_task(WizardTask.RUNNING_SIMULATION):
                self.extra_msg = "Simulation is already in progress, please wait..."
                print(self.extra_msg)
                return

            self.add_task(WizardTask.RUNNING_SIMULATION)
            try:
                if self.sim_type == SimulationType.FULL:
                    self.run_full_simulation()
                else:
                    self.run_partial_simulation()
            except Exception as e:
                print(f"Error while running simulation: {e}.")
            finally:
                self.remove_task(WizardTask.RUNNING_SIMULATION)
                self.extra_msg = ""
                self.update_input_state()

        worker_thread = threading.Thread(
            target=aux,
        )
        worker_thread.start()

    def run_full_simulation(self):
        """Run a full simulation."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        tmp_dir = os.path.join(
            "simulations",
            f"{self.molecule}_{time.strftime('%Y-%m-%d_%H-%M-%S')}_full",
        )
        os.makedirs(tmp_dir)
        cmd.save(os.path.join(tmp_dir, f"{self.molecule}.pdb"), self.molecule)

        simulation = AllAtomSimulationHandler(tmp_dir, self.sim_params)
        simulation.simulate(
            self.molecule,
        )

        output_name = simulation.simulation_data.final_molecule
        cmd.load_traj(os.path.join(tmp_dir, "trajectory.dcd"), output_name)
        cmd.disable(self.molecule)
        cmd.enable(output_name)
        cmd.orient(output_name)

        print(f"Done! Simulation files saved at {tmp_dir}")

    def run_partial_simulation(self):
        """Run a partial simulation."""

        if self.molecule is None:
            print("Please select a molecule.")
            return

        if self.binding_site is None:
            print("Please perform binding site detection first.")
            return

        tmp_dir = os.path.join(
            "simulations",
            f"{self.molecule}_{time.strftime('%Y-%m-%d_%H-%M-%S')}_s{self.sim_params.sim_steps}_r{self.sim_radius}d{self.sim_depth}",
        )
        os.makedirs(tmp_dir)
        cmd.save(os.path.join(tmp_dir, f"{self.molecule}.pdb"), self.molecule)

        simulation = AllAtomSimulationHandler(tmp_dir, self.sim_params)
        simulation.simulate_partial(self.molecule, self.binding_site)

        # constrained_residues = set()
        # for atom in simulation_data.constrained_atoms:
        #     cmd.select(
        #         "constrained_atoms",
        #         f"byres {simulation.simulation_data.final_molecule} and index {atom}",
        #     )
        # if len(simulation_data.constrained_atoms) > 0:
        #     cmd.iterate(
        #         "constrained_atoms",
        #         "constrained_residues.add((resi, chain))",
        #         space=locals(),
        #     )

        output_name = simulation.simulation_data.final_molecule
        cmd.disable(self.molecule)
        cmd.load_traj(os.path.join(tmp_dir, "trajectory.dcd"), output_name)

        self.update_coloring(output_name)
        cmd.align(output_name, self.molecule)
        cmd.disable(self.molecule)
        cmd.enable(output_name)
        cmd.orient(output_name)

        print(f"Done! Simulation files saved at {tmp_dir}")
