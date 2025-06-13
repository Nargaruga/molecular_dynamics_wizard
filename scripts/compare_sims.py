import os
import os.path
from pathlib import Path
import json
import re
import csv

import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from pymol import cmd
import MDAnalysis as mda  # TODO cite the use of rms and align as required by (https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html)
from MDAnalysis.analysis import rms, align

buttons = []

sns.set_theme(style="whitegrid", font_scale=1)


class EmptySelectionError(Exception):
    pass


class ComparisonError(Exception):
    pass


class NoSimulationsError(Exception):
    pass


class RMSF:
    resids: list[int]
    rmsf: list[float]

    def __init__(self, resids, rmsf):
        self.resids = resids
        self.rmsf = rmsf


class Simulation:
    rmsd: list[float]
    rmsd_static: list[float]
    rmsf_by_chain: dict[str, RMSF]
    duration: float
    n_atoms: int
    label: str

    def __init__(self, rmsd, rmsd_static, rmsf_by_chain, duration, n_atoms, label):
        self.rmsd = rmsd
        self.rmsd_static = rmsd_static
        self.rmsf_by_chain = rmsf_by_chain
        self.duration = duration
        self.n_atoms = n_atoms
        self.label = label


class Summary:
    avg_rmsd: float
    avg_rmsd_static: float
    avg_rmsf: float
    duration: float
    n_atoms: int
    radius: int
    depth: int

    def __init__(
        self,
        avg_rmsd,
        avg_rmsd_static,
        avg_rmsf,
        duration,
        n_atoms,
        radius,
        depth,
    ):
        self.avg_rmsd = avg_rmsd
        self.avg_rmsd_static = avg_rmsd_static
        self.avg_rmsf = avg_rmsf
        self.duration = duration
        self.n_atoms = n_atoms
        self.radius = radius
        self.depth = depth


def average_rmsd_over_runs(runs: list[rms.RMSD], n_frames: int) -> list[int]:
    if len(runs) < 1:
        return []

    avg_rmsd = []
    for frame in range(0, n_frames):
        avg_rmsd.append(sum([run.results.rmsd[frame][2] for run in runs]) / len(runs))

    for run in runs:
        print("---RMSD---")
        print(run.results.rmsd)

    return avg_rmsd


def average_rmsd(rmsd: list[float]):
    if not rmsd:
        return []

    return sum(rmsd) / len(rmsd)


def plot_rmsd_static(
    ax: Axes,
    sims: list[Simulation],
    n_frames: int,
):
    for sim in sims:
        if not sim.rmsd_static:
            continue

        ax.plot(
            range(n_frames),
            sim.rmsd_static,
            label=sim.label,
        )

    ax.legend(
        loc="upper left",
        ncols=len(sims),
    )
    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD (Å)")


def get_row(radius):
    if radius == 0:
        row = 0
    elif radius == 4:
        row = 1
    elif radius == 8:
        row = 2
    else:
        raise ComparisonError("Invalid radius")
    return row


def compute_rmsf(topology: str, trajectory: str, selection_str: str) -> RMSF | None:
    if topology == "" or trajectory == "":
        raise ComparisonError("Topology or trajectory file not found")

    if selection_str == "":
        raise EmptySelectionError("No selection string provided")

    u = mda.Universe(topology, trajectory)

    average = align.AverageStructure(
        u, u, select="protein and name CA", ref_frame=0
    ).run()

    ref = average.results.universe
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()

    c_alphas = u.select_atoms(f"protein and name CA and ({selection_str})")

    rmsf = rms.RMSF(c_alphas).run()

    return RMSF(
        list(c_alphas.resids),
        list(rmsf.rmsf),
    )


def average_rmsf_over_runs(
    rmsf_runs_by_chain: dict[str, list[RMSF]],
) -> dict[str, RMSF]:
    if not rmsf_runs_by_chain:
        raise ComparisonError("No runs found")

    avg_rmsf = {}
    for chain, rmsf_runs in rmsf_runs_by_chain.items():
        if not rmsf_runs:
            continue

        avg_rmsf[chain] = RMSF(
            list(rmsf_runs[0].resids), np.zeros_like(rmsf_runs[0].rmsf)
        )
        for rmsf in rmsf_runs:
            avg_rmsf[chain].rmsf += np.array(rmsf.rmsf)

        avg_rmsf[chain].rmsf /= len(rmsf_runs)

    return avg_rmsf


def average_rmsf(rmsf: list[float]) -> float:
    # does this make sense? we are averaging rmsf across all residues

    return sum(rmsf) / len(rmsf)


def annotate_rmsf(ax: Axes, rmsf_by_neigh: list[tuple[list[float], mda.AtomGroup]]):
    if not rmsf_by_neigh:
        return

    for rmsf, c_alphas in rmsf_by_neigh:
        for xy in zip(c_alphas.resids, rmsf):
            ax.annotate("(%s, %.2f)" % xy, xy=xy, textcoords="data")


def toggle_annotation(ax: Axes, rmsf_by_neigh):
    if ax.texts:
        ax.texts.clear()
    else:
        annotate_rmsf(ax, rmsf_by_neigh)

    plt.draw()


def plot_duration(molecule_dir: str):
    molecule_name = Path(molecule_dir).name
    comparison_file = os.path.join(molecule_dir, f"{molecule_name}_comparison.csv")

    data = pd.read_csv(comparison_file)
    data.drop(
        data.columns.difference(["Label", "Duration"]), axis="columns", inplace=True
    )
    data["Label"] = data["Label"].replace("Whole Molecule", "Totale")

    g = sns.barplot(x=data["Label"], y=data["Duration"])

    g.set_xlabel("Tipo di Simulazione")
    g.set_ylabel("Durata (s)")
    for bars in g.containers:
        g.bar_label(bars, fmt="%.1f")

    g.relim()
    g.margins(y=0.08)
    plt.tight_layout()

    plt.savefig(os.path.join(molecule_dir, f"{molecule_name}_duration.png"))
    plt.close("all")


def process_partial_sim(
    molecule_name,
    full_sim,
    minimized,
    dirs,
    radius: int,
    depth: int,
    paratope_selection: str,
    chain_sels: list[str],
    frames: int,
):
    durations = []
    rmsd_runs = []
    rmsd_static_runs = []
    rmsf_runs = {}

    sim_dirs = [f.path for f in os.scandir(dirs) if f.is_dir() and f.name != "ignore"]

    if not sim_dirs:
        return f"No simulation directories found for {molecule_name} with radius {radius} and depth {depth}"

    n_atoms = -1
    for sim_dir in sim_dirs:
        with open(os.path.join(sim_dir, "sim_state.csv")) as f:
            print(f"Processing {sim_dir}")
            final_line = f.readlines()[-1].strip()
            durations.append(float(final_line.split(",")[3]))

        final_molecule_file = os.path.join(
            sim_dir, f"{molecule_name}_minimized_fixed.pdb"
        )

        partial = mda.Universe(
            final_molecule_file,
            os.path.join(sim_dir, "trajectory.dcd"),
        )

        if n_atoms == -1:
            n_atoms = partial.atoms.n_atoms

        rmsd_runs.append(rms.RMSD(partial, full_sim, select=paratope_selection).run())
        rmsd_static_runs.append(
            rms.RMSD(partial, minimized, select=paratope_selection).run(),
        )

        for chain, selection in chain_sels:
            if chain not in rmsf_runs:
                rmsf_runs[chain] = []
            try:
                rmsf_runs[chain].append(
                    compute_rmsf(
                        final_molecule_file,
                        os.path.join(sim_dir, "trajectory.dcd"),
                        selection,
                    )
                )
            except EmptySelectionError:
                continue

    # compute averages over multiple runs
    avg_duration = sum(durations) / max(len(durations), 1)
    avg_rmsd = average_rmsd_over_runs(rmsd_runs, frames)
    avg_rmsd_static = average_rmsd_over_runs(rmsd_static_runs, frames)
    avg_rmsf = average_rmsf_over_runs(rmsf_runs)

    # sim_name = f"sim_r{radius}d{depth}"
    # cmd.load(
    #     os.path.join(sim_dir, f"{molecule_name}_minimized_sliced_fixed.pdb"), sim_name
    # )
    # cmd.load(os.path.join(sim_dir, "trajectory.dcd"), sim_name)

    return (
        avg_duration,
        avg_rmsd,
        avg_rmsd_static,
        avg_rmsf,
        n_atoms,
    )


def get_full_sim_files(molecule_dir, molecule_name, dirs):
    # Find the directory containing the full simulation
    full_sim_regex = re.compile(r".*_full")
    full_sim_dir = None
    for dir in dirs:
        if full_sim_regex.match(dir):
            full_sim_dir = dir
            break

    if not full_sim_dir:
        raise NoSimulationsError("No full simulation directory found")

    # Load the full simulation
    final_full_molecule_file = os.path.join(
        molecule_dir, full_sim_dir, f"{molecule_name}_minimized.pdb"
    )

    final_full_traj_file = os.path.join(molecule_dir, full_sim_dir, "trajectory.dcd")

    sim_data_file = os.path.join(molecule_dir, full_sim_dir, "sim_state.csv")

    return final_full_molecule_file, final_full_traj_file, sim_data_file


def get_paratope_selections(
    paratope_residues, heavy_chains: list[str], light_chains: list[str]
):
    if not paratope_residues:
        raise ComparisonError("No paratope residues found")

    heavy_chains_sels = []
    for hc in heavy_chains:
        heavy_chains_sels.append(
            (
                hc,
                "name CA and ("
                + " or ".join(
                    [
                        f"(resid {res} and chainID {chain})"
                        for res, chain in paratope_residues
                        if chain == hc
                    ]
                )
                + ")",
            ),
        )

    light_chains_sels = []
    for lc in light_chains:
        light_chains_sels.append(
            (
                lc,
                "name CA and ("
                + " or ".join(
                    [
                        f"(resid {res} and chainID {chain})"
                        for res, chain in paratope_residues
                        if chain == lc
                    ]
                )
                + ")",
            )
        )

    paratope = (
        "name CA and ("
        + " or ".join(
            [f"(resid {res} and chainID {chain})" for res, chain in paratope_residues]
        )
        + ")"
    )

    return paratope, heavy_chains_sels, light_chains_sels


def get_sim_duration(sim_data_file):
    with open(sim_data_file) as f:
        final_line = f.readlines()[-1].strip()
        return float(final_line.split(",")[3])


def get_paratope_residues(dirs):
    # quite ugly
    for dir in dirs:
        sim_dirs = [
            f.path for f in os.scandir(dir) if f.is_dir() and f.name != "ignore"
        ]
        for sim_dir in sim_dirs:
            if os.path.exists(os.path.join(sim_dir, "simulated_residues.json")):
                with open(os.path.join(sim_dir, "simulated_residues.json")) as f:
                    sim_metadata = json.load(f)
                    paratope_residues = sim_metadata["paratope"]
                    return paratope_residues

    raise ComparisonError("No paratope residues found in simulation directories")


def serialize_simulations(sims: list[Simulation], base_dir: str, molecule_name: str):
    # Write comparison summary
    with open(
        os.path.join(base_dir, molecule_name, f"{molecule_name}_comparison.csv"), "w"
    ) as f:
        writer = csv.writer(f)
        rmsf_header_labels = []
        chains = sims[0].rmsf_by_chain.keys()
        for chain in chains:
            rmsf_header_labels.append(f"Avg. RMSF {chain}")

        writer.writerow(
            [
                "Label",
                "Duration",
                "Atoms",
                "Avg. RMSD vs Full",
                "Avg. RMSD vs Static",
            ]
            + rmsf_header_labels
        )
        for sim in sims:
            writer.writerow(
                [
                    sim.label,
                    sim.duration,
                    sim.n_atoms,
                    average_rmsd(sim.rmsd),
                    average_rmsd(sim.rmsd_static),
                ]
                + [
                    average_rmsf(sim.rmsf_by_chain[chain].rmsf)
                    for chain in chains
                    if sim.rmsf_by_chain[chain] is not None
                ]
            )

    for sim in sims:
        if sim.label == "Whole Molecule":
            continue

        with open(
            os.path.join(
                base_dir, molecule_name, f"{molecule_name}_{sim.label}_rmsd.csv"
            ),
            "w",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(["Frame", "RMSD"])

            for frame, rmsd in enumerate(sim.rmsd):
                writer.writerow([frame + 1, rmsd])

        with open(
            os.path.join(
                base_dir, molecule_name, f"{molecule_name}_{sim.label}_rmsd_static.csv"
            ),
            "w",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(["Frame", "RMSD"])

            for frame, rmsd in enumerate(sim.rmsd_static):
                writer.writerow([frame + 1, rmsd])

        with open(
            os.path.join(
                base_dir, molecule_name, f"{molecule_name}_{sim.label}_rmsf.csv"
            ),
            "w",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(["chain", "resi", "RMSF"])
            for chain, rmsf in sim.rmsf_by_chain.items():
                for resid, val in zip(rmsf.resids, rmsf.rmsf):
                    writer.writerow([chain, resid, val])


def get_rmsd_plot(molecule_dir: str):
    rmsd_file_regex = re.compile(r".+r([0-9]+)d([0-9]+)_rmsd.csv")

    global_data = pd.DataFrame()
    data_frames = []
    for file in os.listdir(molecule_dir):
        if rmsd_file_regex.match(file):
            radius = int(rmsd_file_regex.match(file).group(1))
            depth = int(rmsd_file_regex.match(file).group(2))
    global_data = pd.DataFrame()
    data_frames = []
    for file in os.listdir(molecule_dir):
        if rmsd_file_regex.match(file):
            radius = int(rmsd_file_regex.match(file).group(1))
            depth = int(rmsd_file_regex.match(file).group(2))

            neigh_data = pd.read_csv(os.path.join(molecule_dir, file))
            downsampled = neigh_data.groupby(neigh_data.index // 5).mean()
            downsampled["Label"] = f"r{radius}d{depth}"
            neigh_data = pd.read_csv(os.path.join(molecule_dir, file))
            downsampled = neigh_data.groupby(neigh_data.index // 5).mean()
            downsampled["Label"] = f"r{radius}d{depth}"

            data_frames.append(downsampled)
            data_frames.append(downsampled)

    global_data = pd.concat(data_frames, ignore_index=True)
    global_data = pd.concat(data_frames, ignore_index=True)

    g = sns.relplot(
        data=global_data,
        kind="line",
        x="Frame",
        y="RMSD",
        hue="Label",
        height=5,
        legend="full",
    )

    molecule_name = Path(molecule_dir).name
    g.ax.set_title(f"{molecule_name}")
    g.ax.set_xlabel("Frame")
    g.ax.set_ylabel("RMSD (Å)")

    return global_data


def get_rmsf_plot(molecule_dir: str):
    molecule_name = Path(molecule_dir).name

    data = pd.read_csv(os.path.join(molecule_dir, f"{molecule_name}_comparison.csv"))
    data = data.drop(
        columns=["Duration", "Atoms", "Avg. RMSD vs Full", "Avg. RMSD vs Static"]
    )
    data = pd.read_csv(os.path.join(molecule_dir, f"{molecule_name}_comparison.csv"))
    data = data.drop(
        columns=["Duration", "Atoms", "Avg. RMSD vs Full", "Avg. RMSD vs Static"]
    )

    # to long format
    data = data.melt(id_vars=["Label"], var_name="Chain", value_name="RMSF")
    data["Chain"] = data["Chain"].str.replace("Avg. RMSF ", "")
    # to long format
    data = data.melt(id_vars=["Label"], var_name="Chain", value_name="RMSF")
    data["Chain"] = data["Chain"].str.replace("Avg. RMSF ", "")

    # we plot the variation in RMSF for the different neighbourhood sizes and each chain
    full_rmsf = data[data["Label"] == "Whole Molecule"].set_index("Chain")["RMSF"]
    data_partial = data[data["Label"] != "Whole Molecule"].copy()
    data_partial["ΔRMSF"] = data_partial.apply(
        lambda row: row["RMSF"] - full_rmsf.get(row["Chain"], 0), axis=1
    )

    heatmap_data = data_partial.pivot(index="Chain", columns="Label", values="ΔRMSF")

    ax = sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".2f",
        center=0.0,
        cmap="coolwarm",
        cbar_kws={"label": "RMSF (Å)"},
    )
    # ax.set_title(f"RMSF of partial simulations for {molecule_name}")
    ax.set_xlabel("Neighbourhood depth")
    ax.set_ylabel("Chain")

    return ax


def plot_size_duration_correlation(base_dir: str):
    # plot a scatterplot of size reduction vs duration gain across all molecules

    global_data = pd.DataFrame()
    data_frames = []
    for molecule_dir in os.listdir(base_dir):
        if molecule_dir == "ignore":
            continue

        molecule_dir = os.path.join(base_dir, molecule_dir)
        if not os.path.isdir(molecule_dir):
            continue

        molecule_name = Path(molecule_dir).name
        data = pd.read_csv(
            os.path.join(molecule_dir, f"{molecule_name}_comparison.csv")
        )

        data["Label"] = data["Label"].replace("Whole Molecule", "Full")
        data.rename(columns={"Label": "Type"}, inplace=True)

        data = data.drop(data.columns.difference(["Type", "Atoms", "Duration"]), axis=1)

        data["Size Reduction"] = (
            (data["Atoms"].iloc[0] - data["Atoms"]) / data["Atoms"].iloc[0]
        ) * 100

        data["Duration Reduction"] = (
            (data["Duration"].iloc[0] - data["Duration"]) / data["Duration"].iloc[0]
        ) * 100

        data_frames.append(data)

    global_data = pd.concat(data_frames, ignore_index=True)

    g = sns.scatterplot(
        data=global_data,
        x="Size Reduction",
        y="Duration Reduction",
        hue="Type",
    )

    # g.set_title("Size reduction vs duration reduction across all molecules")
    g.set_xlabel("Size Reduction (%)")
    g.set_ylabel("Duration Reduction (%)")

    plt.savefig(
        os.path.join(base_dir, "size_vs_duration.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close("all")


def plot_molecule(base_dir: str, molecule_dir: str):
    get_rmsd_plot(molecule_dir)
    plt.savefig(
        os.path.join(base_dir, f"{Path(molecule_dir).name}_rmsd.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.figure()
    get_rmsf_plot(molecule_dir)
    plt.savefig(
        os.path.join(base_dir, f"{Path(molecule_dir).name}_rmsf.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close("all")

    # plt.show()


def plot_all(base_dir: str):
    for molecule_dir in os.listdir(base_dir):
        if molecule_dir == "ignore":
            continue

        molecule_dir = os.path.join(base_dir, molecule_dir)
        if not os.path.isdir(molecule_dir):
            continue

        plot_molecule(base_dir, molecule_dir)

    print("Done!")


def compare_sims(base_dir: str, molecule_name: str, n_frames_str: str):
    n_frames = int(n_frames_str)
    molecule_dir = os.path.join(base_dir, molecule_name)
    chains_file = os.path.join(molecule_dir, "chains_file.json")

    with open(chains_file) as f:
        chains = json.load(f)
        heavy_chains = chains["heavy_chains"]
        light_chains = chains["light_chains"]

    # We have one directory for each simulation
    dirs = [
        f.path for f in os.scandir(molecule_dir) if f.is_dir() and f.name != "ignore"
    ]
    if not dirs:
        print("No simulation directories found")
        return

    final_full_molecule_file, final_full_traj_file, sim_data_file = get_full_sim_files(
        molecule_dir, molecule_name, dirs
    )
    full_sim_system = mda.Universe(
        final_full_molecule_file,
        final_full_traj_file,
    )
    minimized = mda.Universe(final_full_molecule_file)

    paratope, heavy_chain_sels, light_chain_sels = get_paratope_selections(
        get_paratope_residues(dirs), heavy_chains, light_chains
    )

    rmsd_static_full = rms.RMSD(full_sim_system, minimized, select=paratope).run()

    rmsf_full = {}
    for chain, selection in heavy_chain_sels + light_chain_sels:
        rmsf_full[chain] = compute_rmsf(
            final_full_molecule_file,
            final_full_traj_file,
            selection,
        )

    sims = [
        Simulation(
            [],
            list(map(lambda x: x[2], rmsd_static_full.results.rmsd)),
            rmsf_full,
            get_sim_duration(sim_data_file),
            full_sim_system.atoms.n_atoms,
            "Whole Molecule",
        )
    ]

    # cmd.load(final_full_molecule_file, "full_sim")
    # cmd.load(final_full_traj_file, "full_sim")

    for i, dir in enumerate(dirs):
        partial_sim_regex = re.compile("r([0-9]+)d([0-9]+)")
        match = partial_sim_regex.match(Path(dir).name)
        if match:
            radius = match.group(1)
            depth = match.group(2)
        else:
            print(f"Skipping {dir}")
            continue

        (
            avg_duration,
            avg_rmsd,
            avg_rmsd_static,
            avg_rmsf,
            n_atoms,
        ) = process_partial_sim(
            molecule_name,
            full_sim_system,
            minimized,
            dir,
            int(radius),
            int(depth),
            paratope,
            heavy_chain_sels + light_chain_sels,
            n_frames,
        )

        sims.append(
            Simulation(
                avg_rmsd,
                avg_rmsd_static,
                avg_rmsf,
                avg_duration,
                n_atoms,
                f"r{radius}d{depth}",
            )
        )

    serialize_simulations(sims, base_dir, molecule_name)

    print("Done!")


cmd.extend("compare_sims", compare_sims)
cmd.extend("plot_comparison", plot_size_duration_correlation)
cmd.extend("plot_duration", plot_duration)
cmd.extend("plot_molecule", plot_molecule)
cmd.extend("plot_all", plot_all)