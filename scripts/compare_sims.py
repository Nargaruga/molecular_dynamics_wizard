import os
import os.path
from pathlib import Path
import json
import re
import csv

import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import seaborn as sns
import pandas as pd

from pymol import cmd
import MDAnalysis as mda  # TODO cite the use of rms and align as required by (https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html)
from MDAnalysis.analysis import rms, align


sns.set_theme()


buttons = []


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

    return avg_rmsd


def average_rmsd(rmsd: list[float]):
    if not rmsd:
        return []

    return sum(rmsd) / len(rmsd)


# def plot_rmsd(
#     ax: Axes,
#     sims: list[Simulation],
#     n_frames: int,
# ):
#     for sim in sims:
#         if not sim.rmsd:
#             continue

#         ax.plot(
#             range(n_frames),
#             sim.rmsd,
#             label=sim.label,
#         )

#     ax.legend(
#         loc="upper left",
#         ncols=len(sims),
#     )
#     ax.set_xlabel("Frame")
#     ax.set_ylabel("RMSD (Å)")


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


# def serialize_rmsd(rmsd: list[int], n_frames: int, filename: str):
#     if not rmsd or n_frames < 1:
#         return

#     with open(filename, "w") as f:
#         writer = csv.writer(f)
#         writer.writerow(["Frame", "RMSD"])
#         for frame in range(0, n_frames):
#             writer.writerow([frame + 1, rmsd[frame]])


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


def plot_heatmap(rmsf_by_neigh: list[tuple[list[tuple[int, float]], int, int]]):
    """Plot a heatmap of the mean RMSD for a given neighbourhood depth and radius."""

    if not rmsf_by_neigh:
        return

    rmsf_matrix = np.zeros((len(rmsf_by_neigh), len(rmsf_by_neigh[0][0])))
    for rmsf, radius, depth in rmsf_by_neigh:
        row = get_row(radius)
        for resid, val in rmsf:
            np.append(rmsf_matrix[row], val)

    fig, ax = plt.subplots()
    sns.heatmap(
        rmsf_matrix,
        annot=True,
        fmt=".2f",
        ax=ax,
    )

    resids = []
    resids.extend([resid for resid, _ in rmsf_by_neigh[0][0]])

    ax.set_xlabel("Residue ID")
    ax.set_ylabel("Neighbourhood depth")
    ax.set_title("RMSF by neighbourhood depth and radius")
    ax.set_xticklabels(resids)
    ax.set_yticklabels([f"r{radius}d{depth}" for _, radius, depth in rmsf_by_neigh])
    plt.tight_layout()
    plt.show()


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


# def serialize_rmsf(rmsf: tuple[list[float], mda.AtomGroup], filename: str):
#     with open(filename, "w") as f:
#         writer = csv.writer(f)
#         writer.writerow(["resi", "RMSF"])
#         for resi, rmsf_val in zip(rmsf[1].resids, rmsf[0]):
#             writer.writerow([resi, rmsf_val])


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


# def plot_rmsf(ax: Axes, sims: list[Simulation], chain: str):
#     if not sims:
#         return

#     for sim in sims:
#         if chain == "heavy":
#             rmsf = sim.rmsf_h
#         elif chain == "light":
#             rmsf = sim.rmsf_l
#         else:
#             raise ComparisonError("Invalid chain type")

#         if not rmsf:
#             continue

#         ax.plot(
#             rmsf[1].resids,
#             rmsf[0],
#             label=sim.label,
#         )

#     ax.legend(
#         loc="upper left",
#         ncols=len(sims),
#     )
#     ax.set_xlabel("Residue ID")
#     ax.set_ylabel("RMSF (Å)")


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


# def plot_duration(ax: Axes, sims: list[Simulation]):
#     if not sims:
#         return

#     x = np.arange(len(sims))
#     bar_width = 0.05
#     multiplier = 2
#     for sim in sims:
#         offset = bar_width * multiplier
#         rects = ax.bar(x + offset, sim.duration, width=bar_width, label=sim.label)
#         ax.bar_label(rects, fmt="%.2f")
#         multiplier += 1

#     ax.set_ylim(bottom=0)
#     ax.set_xlabel("Neighbourhood depth")
#     ax.set_ylabel("Duration (s)")


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
            final_line = f.readlines()[-1].strip()
            durations.append(float(final_line.split(",")[3]))

        final_molecule_file = os.path.join(
            sim_dir, f"{molecule_name}_minimized_sliced_fixed.pdb"
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


# def parse_simulations_file(file: str) -> list[Summary]:
#     with open(file, "r") as f:
#         reader = csv.reader(f)
#         next(reader)  # Skip the header
#         sims = []
#         for row in reader:
#             label = row[0]
#             if label == "Whole Molecule":
#                 # skip the whole molecule sim for now
#                 continue

#             label_regex = re.compile(r"r([0-9]+)d([0-9]+)")
#             match = label_regex.match(label)
#             if match:
#                 radius = int(match.group(1))
#                 depth = int(match.group(2))
#             else:
#                 # TODO appropriate exception
#                 raise ComparisonError(
#                     f"Invalid label format: {label}. Expected format: r[0-9]+d[0-9]+"
#                 )

#             duration = float(row[1])
#             n_atoms = int(row[2])
#             avg_rmsd = float(row[3])
#             avg_rmsd_static = float(row[4])
#             avg_rmsf_h = float(row[5]) if row[5] else None
#             avg_rmsf_l = float(row[6]) if row[6] else None

#             sim = Summary(
#                 avg_rmsd,
#                 avg_rmsd_static,
#                 avg_rmsf_h,
#                 avg_rmsf_l,
#                 duration,
#                 n_atoms,
#                 radius,
#                 depth,
#             )
#             sims.append(sim)

#     return sims


def plot_rmsd(base_dir: str):
    for molecule_dir in os.listdir(base_dir):
        if molecule_dir == "ignore":
            continue

        molecule_dir = os.path.join(base_dir, molecule_dir)
        if not os.path.isdir(molecule_dir):
            continue

        molecule_name = Path(molecule_dir).name

        rmsd_file_regex = re.compile(r".+r([0-9]+)d([0-9]+)_rmsd.csv")

        global_data = pd.DataFrame()
        data_frames = []
        for file in os.listdir(molecule_dir):
            if rmsd_file_regex.match(file):
                radius = int(rmsd_file_regex.match(file).group(1))
                depth = int(rmsd_file_regex.match(file).group(2))

                neigh_data = pd.read_csv(os.path.join(molecule_dir, file))
                downsampled = neigh_data.groupby(neigh_data.index // 5).mean()
                downsampled["Label"] = f"r{radius}d{depth}"

                data_frames.append(downsampled)

        global_data = pd.concat(data_frames, ignore_index=True)

        g = sns.lineplot(
            data=global_data,
            x="Frame",
            y="RMSD",
            hue="Label",
        )

        g.set_title(f"RMSD of partial simulations for {molecule_name}")
        g.set_xlabel("Frame")
        g.set_ylabel("RMSD (Å)")


def plot_rmsf(base_dir: str):
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
        data = data.drop(
            columns=["Duration", "Atoms", "Avg. RMSD vs Full", "Avg. RMSD vs Static"]
        )

        data = data.melt(
            id_vars=["Label"],
            var_name="Chain",
            value_name="RMSF",
        )

        data["Chain"] = data["Chain"].str.replace("Avg. RMSF ", "")

        g = sns.catplot(
            data=data,
            kind="bar",
            x="Chain",
            y="RMSF",
            hue="Label",
        )

        for container in g.ax.containers:
            g.ax.bar_label(container, fmt="%.2f", padding=2)

        g.ax.set_title(f"RMSF of partial simulations for {molecule_name}")
        g.set_axis_labels("Chain", "RMSF (Å)")


def plot_duration(base_dir: str):
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
        hue="Label",
    )

    g.set_title("Size reduction vs duration reduction across all molecules")
    g.set_xlabel("Size reduction (%)")
    g.set_ylabel("Duration reduction (%)")


def plot_all(dir: str):
    # plot_rmsd(dir)
    # plot_rmsf(dir)
    plot_duration(dir)

    plt.show()


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

    # _, (rmsd_ax, rmsd_static_ax) = plt.subplots(2, 1)
    # rmsd_ax.set_title("Partial VS full simulation")
    # rmsd_static_ax.set_title("Simulation VS static molecule")
    # # compare each frame with the static paratope
    # plot_rmsd_static(rmsd_static_ax, sims, n_frames)
    # # compare each frame of partial simulations with those of the full simulation
    # plot_rmsd(rmsd_ax, sims, n_frames)

    # _, (dur_ax, dur_with_full_ax) = plt.subplots(2, 1)
    # dur_ax.set_title("Simulation duration")
    # # plot durations only for partial simulations since they are much shorter than the full one
    # plot_duration(dur_ax, [sim for sim in sims if sim.label != "Whole Molecule"])
    # # include the full simulation in the runtime plot
    # plot_duration(dur_with_full_ax, sims)

    # fig_rmsf_h, rmsf_ax_h = plt.subplots(layout="constrained")
    # rmsf_ax_h.set_title("RMSF of the paratope on H")
    # # ax_rmsf_h = fig_rmsf_h.add_axes((0.7, 0.05, 0.1, 0.075))
    # # btn_rmsf_h = Button(ax_rmsf_h, "Annotate")
    # # btn_rmsf_h.on_clicked(lambda _: toggle_annotation(rmsf_ax_h, avg_rmsf_h_by_neigh))
    # # buttons.append(btn_rmsf_h)
    # plot_rmsf(rmsf_ax_h, [sim for sim in sims if sim.rmsf_h is not None], "heavy")

    # fig_rmsf_l, rmsf_ax_l = plt.subplots(layout="constrained")
    # rmsf_ax_l.set_title("RMSF of the paratope on L")
    # # ax_rmsf_l = fig_rmsf_l.add_axes((0.7, 0.05, 0.1, 0.075))
    # # btn_rmsf_l = Button(ax_rmsf_l, "Annotate")
    # # btn_rmsf_l.on_clicked(lambda _: toggle_annotation(rmsf_ax_l, avg_rmsf_l_by_neigh))
    # # buttons.append(btn_rmsf_l)
    # plot_rmsf(rmsf_ax_l, [sim for sim in sims if sim.rmsf_l is not None], "light")

    # _, rmsf_ax_full = plt.subplots()
    # rmsf_ax_full.set_title("RMSF on the full molecule")
    # plot_rmsf(
    #     rmsf_ax_full, [sim for sim in sims if sim.label == "Whole Molecule"], "heavy"
    # )

    # plt.show()

    print("Done!")


cmd.extend("compare_sims", compare_sims)
cmd.extend("plot_rmsd", plot_rmsd)
cmd.extend("plot_rmsf", plot_rmsf)
cmd.extend("plot_duration", plot_duration)
cmd.extend("plot_all", plot_all)
