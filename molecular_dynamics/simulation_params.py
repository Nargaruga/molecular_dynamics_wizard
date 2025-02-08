import yaml

from openmm.unit import kelvin, picosecond, femtosecond, nanometer, bar


class SimulationParameters:
    # Preprocessing
    add_solvent: bool = False
    remove_non_simulated: bool = True

    # Minimization
    minimization_steps: int = 1000

    # Equilibration
    nvt_steps: int = 1000
    npt_steps: int = 1000
    eq_pressure: float = 1.0 * bar
    eq_temperature: float = 300.0 * kelvin

    # Simulation
    force_field: str = "charmm36.xml"
    water_model: str = "charmm36/water.xml"
    sim_steps: int = 100000

    # Integrator
    temperature: float = 300.0 * kelvin
    friction_coeff: float = 1.0 / picosecond
    timestep: float = 2.0 * femtosecond

    # Reporters
    report_interval: int = 1000

    def print(self):
        print(f"Add solvent: {self.add_solvent}")
        print(f"Remove non-simulated: {self.remove_non_simulated}")
        print(f"Minimization steps: {self.minimization_steps}")
        print(f"NVT steps: {self.nvt_steps}")
        print(f"NPT steps: {self.npt_steps}")
        print(f"Equilibration pressure: {self.eq_pressure}")
        print(f"Equilibration temperature: {self.eq_temperature}")
        print(f"Force field: {self.force_field}")
        print(f"Water model: {self.water_model}")
        print(f"Simulation steps: {self.sim_steps}")
        print(f"Temperature: {self.temperature}")
        print(f"Friction coefficient: {self.friction_coeff}")
        print(f"Timestep: {self.timestep}")
        print(f"Report interval: {self.report_interval}")

    def parse_file(self, file_name):
        with open(file_name, "r") as file:
            try:
                config = yaml.safe_load(file)

                self.add_solvent = config["add_solvent"]
                self.remove_non_simulated = config["remove_non_simulated"]

                self.minimization_steps = config["minimization_steps"]

                self.nvt_steps = config["nvt_steps"]
                self.npt_steps = config["npt_steps"]
                self.eq_pressure = config["eq_pressure"] * bar
                self.eq_temperature = config["eq_temperature"] * kelvin

                self.force_field = config["force_field"]
                self.water_model = config["water_model"]
                self.sim_steps = config["sim_steps"]

                self.temperature = config["temperature"] * kelvin
                self.friction_coeff = config["friction_coeff"] / picosecond
                self.timestep = config["timestep"] * femtosecond

                self.report_interval = config["report_interval"]
            except yaml.YAMLError:
                print(f'Error reading configuration file "{file_name}"')
