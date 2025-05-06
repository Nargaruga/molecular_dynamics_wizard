from pymol import cmd, CmdException


class BindingSiteError(Exception):
    """Custom exception for binding site errors."""

    pass


def select_neighbourhood(new_sel, base_sel, allowed_sel, radius=0, depth=1):
    # Expand the base selection by a given radius
    try:
        if radius > 0:
            cmd.select(
                name=new_sel,
                selection=f"byres {base_sel} expand {radius} and ({allowed_sel})",
                merge=1,
            )
        else:
            cmd.select(name=new_sel, selection=base_sel)

        # Extend the newly created selection by a given depth
        if depth > 0:
            cmd.select(
                name=new_sel,
                selection=f"byres {new_sel} extend {depth} and ({allowed_sel})",
                merge=1,
            )
    except CmdException as e:
        raise BindingSiteError(f"Error selecting neighbourhood: {e}") from e


def remove_overlap(sel1, sel2):
    try:
        cmd.select(
            name=sel1,
            selection=f"byres {sel1} and not byres {sel2}",
        )
    except CmdException as e:
        raise BindingSiteError(f"Error removing overlap: {e}") from e


def get_residues(selection_name):
    try:
        residues = set()
        cmd.iterate(
            f"{selection_name}",
            "residues.add((resi, chain))",
            space={"residues": residues},
        )
        return residues

    except CmdException as e:
        raise BindingSiteError(f"Error getting residues: {e}") from e


def residues_to_atoms(molecule: str, residues: set) -> set:
    try:
        atoms = set()
        for residue in residues:
            cmd.iterate(
                f"{molecule} and resi {residue[0]} and chain {residue[1]}",
                "atoms.add(index)",
                space={"atoms": atoms},
            )
        return atoms

    except CmdException as e:
        raise BindingSiteError(f"Error converting residues to atoms: {e}") from e


class BindingSite:
    def __init__(self, molecule, heavy_chains, light_chains):
        self.molecule = molecule
        self.heavy_chains = heavy_chains
        self.light_chains = light_chains
        self.paratope_sel = f"{molecule}_par"
        self.epitope_sel = f"{molecule}_epi"
        self.paratope_neigh_sel = f"{molecule}_par_neigh"
        self.epitope_neigh_sel = f"{molecule}_epi_neigh"
        self.ext_paratope_neigh_sel = f"{molecule}_ext_par_neigh"
        self.ext_epitope_neigh_sel = f"{molecule}_ext_epi_neigh"

    def reset(self):
        """Reset the binding site selections."""
        try:
            cmd.delete(self.paratope_sel)
            cmd.delete(self.epitope_sel)
            cmd.delete(self.paratope_neigh_sel)
            cmd.delete(self.epitope_neigh_sel)
            cmd.delete(self.ext_paratope_neigh_sel)
            cmd.delete(self.ext_epitope_neigh_sel)
        except CmdException as e:
            raise BindingSiteError(f"Error deleting selections: {e}") from e

    def update_neighbourhoods(self, radius, depth):
        try:
            cmd.delete(self.paratope_neigh_sel)
            cmd.delete(self.epitope_neigh_sel)
            cmd.delete(self.ext_paratope_neigh_sel)
            cmd.delete(self.ext_epitope_neigh_sel)
        except CmdException as e:
            raise BindingSiteError(f"Error deleting selections: {e}") from e

        self.select_paratope_neigh(radius, depth)
        self.select_epitope_neigh(radius, depth)
        self.select_ext_paratope_neigh()
        self.select_ext_epitope_neigh()

        remove_overlap(self.ext_paratope_neigh_sel, self.paratope_neigh_sel)
        remove_overlap(self.paratope_neigh_sel, self.paratope_sel)

        remove_overlap(self.ext_epitope_neigh_sel, self.epitope_neigh_sel)
        remove_overlap(self.epitope_neigh_sel, self.epitope_sel)

    def select_paratope(self):
        """Identify the paratope on the selected antibody through the appropriate wizard."""

        print(
            f"Selecting paratope for {self.molecule}, chains {self.heavy_chains} and {self.light_chains}..."
        )
        try:
            cmd.wizard("paratope")
            cmd.get_wizard().set_molecule(self.molecule)

            for heavy_chain in self.heavy_chains:
                if not cmd.get_wizard().is_heavy_chain_selected(heavy_chain):
                    cmd.get_wizard().toggle_heavy_chain(heavy_chain)
            for light_chain in self.light_chains:
                if not cmd.get_wizard().is_light_chain_selected(light_chain):
                    cmd.get_wizard().toggle_light_chain(light_chain)

            cmd.get_wizard().set_selection_name(self.paratope_sel)
            cmd.get_wizard().set_highlight(False)
            cmd.get_wizard().run(block=True)
            cmd.get_wizard().toggle_label_pos()
            cmd.set_wizard()
        except CmdException as e:
            raise BindingSiteError(f"Error selecting paratope: {e}") from e

    def select_epitope(self):
        """Get the residues of the epitope based on the paratope selection."""

        try:
            cmd.select(
                name=self.epitope_sel,
                selection=f"byres {self.molecule} and not ("
                + (
                    " or ".join(
                        [
                            f"chain {chain}"
                            for chain in self.heavy_chains + self.light_chains
                        ]
                    )
                )
                + f") near_to 6.0 of {self.paratope_sel}",
            )
        except CmdException as e:
            raise BindingSiteError(f"Error selecting epitope: {e}") from e

    def select_paratope_neigh(self, radius, depth):
        allowed_sel = (
            f"byres {self.molecule} and ("
            + (
                " or ".join(
                    [
                        f"chain {chain}"
                        for chain in self.heavy_chains + self.light_chains
                    ]
                )
            )
            + ")"
        )

        select_neighbourhood(
            self.paratope_neigh_sel,
            self.paratope_sel,
            allowed_sel,
            radius,
            depth,
        )

    def select_epitope_neigh(
        self,
        radius,
        depth,
    ):
        allowed_sel = (
            f"byres {self.molecule} and ("
            + (
                " and ".join(
                    [
                        f"not chain {chain}"
                        for chain in self.heavy_chains + self.light_chains
                    ]
                )
            )
            + ")"
        )

        select_neighbourhood(
            self.epitope_neigh_sel,
            self.epitope_sel,
            allowed_sel,
            radius,
            depth,
        )

    def select_ext_paratope_neigh(self):
        allowed_sel = (
            f"byres {self.molecule} and ("
            + (
                " or ".join(
                    [
                        f"chain {chain}"
                        for chain in self.heavy_chains + self.light_chains
                    ]
                )
            )
            + ")"
        )

        select_neighbourhood(
            self.ext_paratope_neigh_sel,
            self.paratope_neigh_sel,
            allowed_sel,
        )

    def select_ext_epitope_neigh(self):
        allowed_sel = (
            f"byres {self.molecule} and ("
            + (
                " and ".join(
                    [
                        f"not chain {chain}"
                        for chain in self.heavy_chains + self.light_chains
                    ]
                )
            )
            + ")"
        )

        select_neighbourhood(
            self.ext_epitope_neigh_sel,
            self.epitope_neigh_sel,
            allowed_sel,
        )
