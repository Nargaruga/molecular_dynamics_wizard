import sys
import os
import json

from openmm.app.pdbfile import PDBFile
from openmm.app import Simulation
from openmm.openmm import (
    LangevinIntegrator,
    MonteCarloBarostat,
)
from openmm.app import ForceField, CutoffNonPeriodic, StateDataReporter, DCDReporter
from openmm.unit import (
    kelvin,
    picosecond,
    femtosecond,
    nanometer,
    bar,
)

import pdbfixer

from pymol import cmd

from .simulation_params import SimulationParameters
from .binding_site import BindingSite, get_residues, residues_to_atoms


class SimulationData:
    def __init__(self):
        self.final_molecule = ""
        self.constrained_atoms = set()
        self.binding_site = None


class AllAtomSimulationHandler:
    def __init__(self, tmp_dir, parameters: SimulationParameters):
        self.tmp_dir = tmp_dir
        self.parameters = parameters
        self.simulation_data = SimulationData()

        print("Loaded simulation parameters:")
        parameters.print()

    def fix_pdb(self, input_name, output_name):
        fixer = pdbfixer.PDBFixer(
            filename=os.path.join(self.tmp_dir, f"{input_name}.pdb")
        )
        fixer.findMissingResidues()
        chains = list(fixer.topology.chains())
        keys_to_remove = []
        for key in fixer.missingResidues.keys():
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                keys_to_remove.append(key)
        for key in keys_to_remove:
            del fixer.missingResidues[key]
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

    def slice_object(self, molecule: str, residues_to_keep: set, output_name: str):
        selection = ""
        for resi, chain in residues_to_keep:
            selection = "selection_to_keep"
            cmd.select(
                "selection_to_keep",
                f"{molecule} and resi {resi} and chain {chain}",
                merge=1,
            )

        cmd.create(output_name, selection)
        cmd.delete(selection)

    def create_system(self, molecule):
        pdb = PDBFile(os.path.join(self.tmp_dir, f"{molecule}.pdb"))
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

        return pdb, system, integrator

    def minimize(self, molecule, output_name):
        pdb, system, integrator = self.create_system(molecule)

        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=self.parameters.minimization_steps)

        self.snapshot(simulation, f"{output_name}.pdb")

    def simulate(self, molecule):
        print("Fixing PDB...")
        final_molecule = f"{molecule}_fixed"
        self.fix_pdb(molecule, final_molecule)
        self.simulation_data.final_molecule = final_molecule

        pdb, system, integrator = self.create_system(final_molecule)
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

    def simulate_partial(self, molecule, binding_site: BindingSite):
        residues_to_simulate = (
            get_residues(binding_site.paratope_sel)
            .union(get_residues(binding_site.paratope_neigh_sel))
            .union(get_residues(binding_site.epitope_sel))
            .union(get_residues(binding_site.epitope_neigh_sel))
        )

        # Create a new object out of the non-simulated atoms
        non_sim_molecule = f"{molecule}_non_sim"
        non_sim_molecule_sel = f"{non_sim_molecule}_sel"
        # This object includes the locked neighbourhood, which will help us
        # align the final trajaectory with the original structure
        cmd.select(
            name=non_sim_molecule_sel,
            selection=f"byres {molecule}"
            + f" and not {binding_site.paratope_sel}"
            + f" and not {binding_site.paratope_neigh_sel}"
            + f" and not {binding_site.epitope_sel}"
            + f" and not {binding_site.epitope_neigh_sel}",
        )
        cmd.create(non_sim_molecule, non_sim_molecule_sel)
        cmd.delete(non_sim_molecule_sel)
        cmd.align(non_sim_molecule, molecule)

        if self.parameters.remove_non_simulated:
            print("Removing non-simulated atoms...")
            # This operation causes the renumbering of all atoms,
            # but the residue ids are preserved
            residues_to_keep = residues_to_simulate.union(
                get_residues(binding_site.ext_paratope_neigh_sel)
            ).union(get_residues(binding_site.ext_epitope_neigh_sel))

            sliced_molecule = f"{molecule}_sliced"
            self.slice_object(
                molecule,
                residues_to_keep,
                sliced_molecule,
            )

            to_fix = sliced_molecule
        else:
            to_fix = molecule

        cmd.save(os.path.join(self.tmp_dir, f"{to_fix}.pdb"), to_fix)
        cmd.disable(to_fix)

        print("Fixing PDB...")
        fixed_molecule = f"{to_fix}_fixed"
        self.fix_pdb(to_fix, fixed_molecule)
        cmd.load(os.path.join(self.tmp_dir, f"{fixed_molecule}.pdb"), fixed_molecule)
        cmd.disable(fixed_molecule)

        atoms_to_simulate = residues_to_atoms(
            fixed_molecule,
            residues_to_simulate,
        )

        # Save residue ids to json
        with open(os.path.join(self.tmp_dir, "simulated_residues.json"), "w") as f:
            json.dump(
                {
                    "paratope": list(get_residues(binding_site.paratope_sel)),
                    "paratope_neighbourhood": list(
                        get_residues(binding_site.paratope_neigh_sel)
                    ),
                    "locked_paratope_neighbourhood": list(
                        get_residues(binding_site.ext_paratope_neigh_sel)
                    ),
                    "epitope": list(get_residues(binding_site.epitope_sel)),
                    "epitope_neighbourhood": list(
                        get_residues(binding_site.epitope_neigh_sel)
                    ),
                    "locked_epitope_neighbourhood": list(
                        get_residues(binding_site.ext_epitope_neigh_sel)
                    ),
                },
                f,
            )

        pdb, system, integrator = self.create_system(fixed_molecule)

        for i in range(system.getNumConstraints()):
            particle1, particle2, _ = system.getConstraintParameters(i)
            self.simulation_data.constrained_atoms.add(particle1)
            self.simulation_data.constrained_atoms.add(particle2)

        for atom in pdb.topology.atoms():
            if (
                atom.index not in atoms_to_simulate
                and atom.index not in self.simulation_data.constrained_atoms
            ):
                system.setParticleMass(atom.index, 0.0)

        simulation = Simulation(pdb.topology, system, integrator)
        self.enable_reporters(simulation)
        simulation.context.setPositions(pdb.positions)

        print("Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=self.parameters.minimization_steps)
        # Explicitly write the conect records. We will need them when we
        # load the minimized molecule in connect_mode 1
        minimized_molecule = f"{fixed_molecule}_minimized"
        # cmd.set("pdb_conect_all", 1)
        self.snapshot(simulation, f"{minimized_molecule}.pdb")
        # cmd.set("pdb_conect_all", 0)

        self.simulation_data.final_molecule = minimized_molecule
        self.simulation_data.binding_site = binding_site

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
