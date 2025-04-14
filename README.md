# Molecular Dynamics Wizard
PyMOL wizard for performing molecular dynamics simulations with OpenMM. Features the ability to limit the scope of a simulation to the binding site of an antibody-antigen complex.

## Installation
The wizard can be installed with the [Wizard Installer](https://github.com/Nargaruga/pymol_wizard_installer) tool.
The installation process will also create a `plugin.zip` file which can be used to install the configuration plugin through PyMOL's plugin manager.

## Usage
The wizard can be accessed in one of two ways:
- by writing `wizard dynamics` in the PyMOL console;
- by going to `Wizard->Molecular Dynamics` in the external GUI (or the internal one, if in VR mode);

### Minimizing a Structure
You can minimize the molecule selected through the drop-down menu with the `Minimize Energy` button.

### Running Simulations
After selecting a molecule through the menu, you are given two choices for the simulation type:
- `Full`: perform a molecular dynamics simulation of the whole molecule;
- `Partial`: restrict the simulation to the binding site and its surrounding area. Only available for antibodies or antibody-antigen complexes.

Running a partial simulation requires you to select the antibody's heavy and light chain and use `Detect Binding Site`. This will find the paratope and, if applicable, the epitope, which is defined as any residue not belonging to the heavy or light chain found in a 6Å radius around the paratope. Once that's done, the molecule is colored as follows:
- green (yellow): paratope (epitope) residues;
- blue (cyan): residues neighbouring the paratope (epitope);
- red (magenta): residues bordering the paratope (epitope) neighbourhood.
Any residue not included in the above is removed from the molecule. Residues in the red (magenta) set are locked in place during the simulation. You can adjust the neighbourhood radius and depth as required to increase or decrease the portion of molecule to be simulated. The `Run Simulation` button starts the simulation in the background, automatically loading the trajectory as soon as it is ready.

### Configuration
Simulation parameters can be configured through the companion plugin.

