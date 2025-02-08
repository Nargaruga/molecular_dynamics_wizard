import os
import sys
from abc import ABC, abstractmethod

from openmm.app.pdbxfile import PDBFile
from openmm.app import StateDataReporter, DCDReporter


class SimulationHandler(ABC):
    def __init__(self, tmp_dir, parameters):
        self.tmp_dir = tmp_dir
        self.parameters = parameters

    @abstractmethod
    def preprocess_input(self, input_molecule, output_name):
        """Prepare the input files for the simulation."""
        pass

    @abstractmethod
    def simulate(self, molecule, atoms_to_simulate):
        """Perform a molecular dynamics simulation limited to the specified atoms."""
        pass

    @abstractmethod
    def postprocess_output(self):
        pass

    def snapshot(self, simulation, file_name):
        """Save a snapshot of the simulation."""

        state = simulation.context.getState(getPositions=True)
        with open(os.path.join(self.tmp_dir, file_name), "w") as output:
            PDBFile.writeFile(
                simulation.topology, state.getPositions(), output, keepIds=True
            )

    def enable_reporters(self, simulation):
        """Enable the reporters for the simulation."""

        simulation.reporters.append(
            DCDReporter(
                os.path.join(self.tmp_dir, "trajectory.dcd"),
                self.parameters.report_interval,
            )
        )

        simulation.reporters.append(
            StateDataReporter(
                os.path.join(self.tmp_dir, "sim_state.csv"),
                self.parameters.report_interval,
                step=True,
                elapsedTime=True,
                potentialEnergy=True,
                temperature=True,
            )
        )

        simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                self.parameters.report_interval,
                step=True,
                totalSteps=self.parameters.sim_steps,
                elapsedTime=True,
                progress=True,
                remainingTime=True,
            )
        )
