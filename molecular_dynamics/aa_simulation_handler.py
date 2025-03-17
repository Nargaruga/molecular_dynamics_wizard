import os

from openmm.app.pdbfile import PDBFile
from openmm.app import Simulation
from openmm.openmm import LangevinIntegrator, MonteCarloBarostat
from openmm.app import ForceField, CutoffNonPeriodic, HBonds
from openmm.unit import kelvin, picosecond, femtosecond, nanometer, bar, amu

import pdbfixer

from pymol import cmd

from .simulation_handler import SimulationHandler
from .simulation_params import SimulationParameters


class AllAtomSimulationHandler(SimulationHandler):
    def __init__(self, tmp_dir, parameters: SimulationParameters):
        self.tmp_dir = tmp_dir
        self.parameters = parameters

        print("Loaded simulation parameters:")
        parameters.print()

    def preprocess_input(self, input_molecule, output_name):
        """Prepare the input files for the simulation."""

        fixer = pdbfixer.PDBFixer(
            filename=os.path.join(self.tmp_dir, f"{input_molecule}.pdb")
        )
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens()
        if self.parameters.add_solvent:
            fixer.addSolvent(padding=1.0 * nanometer)

        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            open(os.path.join(self.tmp_dir, output_name), "w"),
            keepIds=True,
        )

    def residues_to_atoms(self, molecule: str, residues: set) -> set:
        atoms = set()
        for residue in residues:
            cmd.iterate(
                f"{molecule} and resi {residue[0]} and chain {residue[1]}",
                "atoms.add(index)",
                space=locals(),
            )
        return atoms

    def slice_object(self, molecule: str, residues_to_keep: set):
        selection = ""
        for resi, chain in residues_to_keep:
            selection = "selection_to_keep"
            cmd.select(
                "selection_to_keep",
                f"{molecule} and resi {resi} and chain {chain}",
                merge=1,
            )

        cmd.create(f"{molecule}_sliced", selection)
        cmd.delete(selection)

    def simulate(self, molecule, atoms_to_simulate=None):
        """Perform a molecular dynamics simulation limited to the specified atoms."""

        pdb = PDBFile(os.path.join(self.tmp_dir, f"{molecule}.pdb"))
        forcefield = ForceField(
            self.parameters.force_field, self.parameters.water_model
        )

        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=CutoffNonPeriodic,
            nonbondedCutoff=2.0 * nanometer,
            constraints=HBonds,
            hydrogenMass=1.5 * amu,
        )

        integrator = LangevinIntegrator(
            self.parameters.temperature * kelvin,
            self.parameters.friction_coeff / picosecond,
            self.parameters.timestep * femtosecond,
        )

        constrained_atoms = set()
        for i in range(system.getNumConstraints()):
            particle1, particle2, _ = system.getConstraintParameters(i)
            constrained_atoms.add(particle1)
            constrained_atoms.add(particle2)

        if atoms_to_simulate is not None:
            for atom in pdb.topology.atoms():
                if (
                    atom.index not in atoms_to_simulate
                    and atom.index not in constrained_atoms
                ):
                    system.setParticleMass(atom.index, 0)

        simulation = Simulation(pdb.topology, system, integrator)
        self.enable_reporters(simulation)
        simulation.context.setPositions(pdb.positions)

        print("Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=self.parameters.minimization_steps)

        self.snapshot(simulation, f"{molecule}_minimized.pdb")

        if self.parameters.nvt_steps > 0:
            print("NVT Equilibration...")
            simulation.step(self.parameters.nvt_steps)

        if self.parameters.npt_steps > 0:
            print("NPT Equilibration...")
            system.addForce(
                MonteCarloBarostat(
                    self.parameters.eq_pressure * bar,
                    self.parameters.eq_temperature * kelvin,
                )
            )
            simulation.context.reinitialize(preserveState=True)
            simulation.step(self.parameters.npt_steps)

        print("Running simulation...")
        simulation.step(self.parameters.sim_steps)

        return constrained_atoms

    def postprocess_output(self):
        pass
