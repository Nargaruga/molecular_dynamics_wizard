import sys
import os
import json

from openmm.app.pdbfile import PDBFile
from openmm.app import Simulation
from openmm.openmm import (
    LangevinIntegrator,
    MonteCarloBarostat,
    CustomExternalForce,
)
from openmm.app import ForceField, CutoffNonPeriodic, StateDataReporter, DCDReporter
from openmm.unit import (
    kelvin,
    picosecond,
    femtosecond,
    nanometer,
    bar,
    kilojoules_per_mole,
)

import pdbfixer

from pymol import cmd

from .simulation_params import SimulationParameters


class SimulationData:
    def __init__(self):
        self.final_molecule = ""
        self.constrained_atoms = set()

        self.paratope_residues = set()
        self.paratope_neigh_residues = set()
        self.locked_paratope_neigh_residues = set()

        self.epitope_residues = set()
        self.epitope_neigh_residues = set()
        self.locked_epitope_neigh_residues = set()


class AllAtomSimulationHandler:
    def __init__(self, tmp_dir, parameters: SimulationParameters):
        self.tmp_dir = tmp_dir
        self.parameters = parameters
        self.simulation_data = SimulationData()

        print("Loaded simulation parameters:")
        parameters.print()

    def shrink_to_interaction_site(self, molecule, depth, heavy_chains, light_chains):
        if depth >= 0:
            self.identify_paratope(
                molecule, f"{molecule}_paratope", heavy_chains, light_chains
            )

            self.select_epitope(
                molecule,
                f"{molecule}_paratope",
                f"{molecule}_epitope",
                heavy_chains,
                light_chains,
            )

            cmd.iterate(
                f"{molecule}_paratope",
                "self.simulation_data.paratope_residues.add((resi, chain))",
                space=locals(),
            )

            cmd.iterate(
                f"{molecule}_epitope",
                "self.simulation_data.epitope_residues.add((resi, chain))",
                space=locals(),
            )

            # Get the neighbourhood of the paratope
            self.simulation_data.paratope_neigh_residues = self.get_neigbourhood(
                f"{molecule}_paratope_neigh",
                f"{molecule}_paratope",
                " or ".join(
                    [f"chain {chain}" for chain in heavy_chains + light_chains]
                ),
                depth,
            )

            # Obtain a further extended neighbourhood of the paratope that will have its atoms locked
            self.simulation_data.locked_paratope_neigh_residues = self.get_neigbourhood(
                f"{molecule}_locked_paratope_neigh",
                f"{molecule}_paratope_neigh",
                " or ".join(
                    [f"chain {chain}" for chain in heavy_chains + light_chains]
                ),
            )

            # Get the neighbourhood of the epitope
            self.simulation_data.epitope_neigh_residues = self.get_neigbourhood(
                f"{molecule}_epitope_neigh",
                f"{molecule}_epitope",
                " and ".join(
                    [f"not chain {chain}" for chain in heavy_chains + light_chains]
                ),
                depth,
            )

            # Obtain a further extended neighbourhood of the epitope that will have its atoms locked
            self.simulation_data.locked_epitope_neigh_residues = self.get_neigbourhood(
                f"{molecule}_locked_epitope_neigh",
                f"{molecule}_epitope_neigh",
                " and ".join(
                    [f"not chain {chain}" for chain in heavy_chains + light_chains]
                ),
            )

            # Create a new object out of the non-simulated atoms
            non_sim_molecule = f"{molecule}_non_sim"
            non_sim_molecule_sel = f"{non_sim_molecule}_sel"
            # This object includes the locked neighbourhood, which will help us
            # align the final trajaectory with the original structure
            cmd.select(
                name=non_sim_molecule_sel,
                selection=f"byres {molecule} and not {molecule}_paratope_neigh and not {molecule}_epitope_neigh",
            )
            cmd.create(non_sim_molecule, non_sim_molecule_sel)
            cmd.delete(non_sim_molecule_sel)
            cmd.align(non_sim_molecule, molecule)
            cmd.disable(non_sim_molecule)

            if self.parameters.remove_non_simulated:
                print("Removing non-simulated atoms...")
                # This operation causes the renumbering of all atoms,
                # but the residue ids are preserved
                self.slice_object(
                    molecule,
                    self.simulation_data.locked_paratope_neigh_residues.union(
                        self.simulation_data.locked_epitope_neigh_residues
                    ),
                )

                to_fix = f"{molecule}_sliced"
            else:
                to_fix = molecule
        else:
            to_fix = molecule

        cmd.save(os.path.join(self.tmp_dir, f"{to_fix}.pdb"), to_fix)
        cmd.disable(to_fix)

        final_molecule = f"{molecule}_fixed"
        self.fix_pdb(to_fix, final_molecule)

        return final_molecule

    def fix_pdb(self, input_name, output_name):
        fixer = pdbfixer.PDBFixer(
            filename=os.path.join(self.tmp_dir, f"{input_name}.pdb")
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
            open(os.path.join(self.tmp_dir, f"{output_name}.pdb"), "w"),
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

    def get_neigbourhood(self, selection_name, target_name, chains, depth=1):
        residues = set()

        if depth == 0:
            cmd.select(name=selection_name, selection=target_name)
        else:
            cmd.select(
                name=selection_name,
                selection=f"byres {target_name} extend {depth} and ({chains})",
                merge=1,
            )

        cmd.iterate(
            selection_name,
            "residues.add((resi, chain))",
            space=locals(),
        )

        return residues

    def identify_paratope(self, molecule, selection_name, heavy_chains, light_chains):
        """Identify the paratope on the selected antibody through the appropriate wizard."""

        cmd.wizard("paratope")
        cmd.get_wizard().set_molecule(molecule)
        for heavy_chain in heavy_chains:
            cmd.get_wizard().set_heavy_chain(heavy_chain)
        for light_chain in light_chains:
            cmd.get_wizard().set_light_chain(light_chain)
        cmd.get_wizard().set_selection_name(selection_name)
        cmd.get_wizard().run()
        cmd.get_wizard().toggle_label_pos()
        cmd.set_wizard()

    def select_epitope(
        self, molecule, paratope_sel, selection_name, heavy_chains, light_chains
    ):
        """Get the residues of the epitope based on the paratope selection."""

        epitope_sel = selection_name
        cmd.select(
            name=epitope_sel,
            selection=f"byres {molecule} and not ("
            + (" or ".join([f"chain {chain}" for chain in heavy_chains + light_chains]))
            + f") near_to 6.0 of {paratope_sel}",
        )

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

    def simulate(self, molecule):
        final_molecule = f"{molecule}_fixed"
        self.fix_pdb(molecule, final_molecule)

        pdb = PDBFile(os.path.join(self.tmp_dir, f"{final_molecule}.pdb"))
        forcefield = ForceField(
            self.parameters.force_field, self.parameters.water_model
        )

        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=CutoffNonPeriodic,
            nonbondedCutoff=2.0 * nanometer,
        )

        integrator = LangevinIntegrator(
            self.parameters.temperature * kelvin,
            self.parameters.friction_coeff / picosecond,
            self.parameters.timestep * femtosecond,
        )

        for i in range(system.getNumConstraints()):
            particle1, particle2, _ = system.getConstraintParameters(i)
            self.simulation_data.constrained_atoms.add(particle1)
            self.simulation_data.constrained_atoms.add(particle2)

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

    def simulate_partial(self, molecule, depth, heavy_chains, light_chains):
        if depth < 0:
            print("Neighbourhood depth must be >= 0.")
            return

        final_molecule = self.shrink_to_interaction_site(
            molecule, depth, heavy_chains, light_chains
        )
        self.simulation_data.final_molecule = final_molecule
        cmd.disable(molecule)
        cmd.load(os.path.join(self.tmp_dir, f"{final_molecule}.pdb"))
        cmd.show_as("licorice", final_molecule)

        # Remove overlap
        self.simulation_data.locked_paratope_neigh_residues = (
            self.simulation_data.locked_paratope_neigh_residues
            - self.simulation_data.paratope_neigh_residues
        )
        self.simulation_data.paratope_neigh_residues = (
            self.simulation_data.paratope_neigh_residues
            - self.simulation_data.paratope_residues
        )

        self.simulation_data.locked_epitope_neigh_residues = (
            self.simulation_data.locked_epitope_neigh_residues
            - self.simulation_data.epitope_neigh_residues
        )
        self.simulation_data.epitope_neigh_residues = (
            self.simulation_data.epitope_neigh_residues
            - self.simulation_data.epitope_residues
        )

        # Save residue ids to json
        with open(os.path.join(self.tmp_dir, "simulated_residues.json"), "w") as f:
            json.dump(
                {
                    "paratope": list(self.simulation_data.paratope_residues),
                    "paratope_neighbourhood": list(
                        self.simulation_data.paratope_neigh_residues
                    ),
                    "locked_paratope_neighbourhood": list(
                        self.simulation_data.locked_paratope_neigh_residues
                    ),
                    "epitope": list(self.simulation_data.epitope_residues),
                    "epitope_neighbourhood": list(
                        self.simulation_data.epitope_neigh_residues
                    ),
                    "locked_epitope_neighbourhood": list(
                        self.simulation_data.locked_epitope_neigh_residues
                    ),
                    "depth": depth,
                },
                f,
            )

        pdb = PDBFile(os.path.join(self.tmp_dir, f"{final_molecule}.pdb"))
        forcefield = ForceField(
            self.parameters.force_field, self.parameters.water_model
        )

        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=CutoffNonPeriodic,
            nonbondedCutoff=2.0 * nanometer,
        )

        integrator = LangevinIntegrator(
            self.parameters.temperature * kelvin,
            self.parameters.friction_coeff / picosecond,
            self.parameters.timestep * femtosecond,
        )

        for i in range(system.getNumConstraints()):
            particle1, particle2, _ = system.getConstraintParameters(i)
            self.simulation_data.constrained_atoms.add(particle1)
            self.simulation_data.constrained_atoms.add(particle2)

        atoms_to_simulate = self.residues_to_atoms(
            final_molecule,
            (
                self.simulation_data.paratope_residues.union(
                    self.simulation_data.paratope_neigh_residues
                )
                .union(self.simulation_data.epitope_residues)
                .union(self.simulation_data.epitope_neigh_residues)
            ),
        )

        restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        system.addForce(restraint)
        restraint.addGlobalParameter("k", 1000.0 * kilojoules_per_mole / nanometer)
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")

        for atom in pdb.topology.atoms():
            if (
                atom.index not in atoms_to_simulate
                and atom.index not in self.simulation_data.constrained_atoms
            ):
                restraint.addParticle(atom.index, pdb.positions[atom.index])

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
