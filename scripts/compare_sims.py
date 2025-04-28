import os
import os.path
import json
import re
import csv

import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

from pymol import cmd
import MDAnalysis as mda  # TODO cite the use of rms and align as required by (https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html)
from MDAnalysis.analysis import rms, align

buttons = []


class NoSimulationsException(Exception):
    def __init__(self, msg):
        self.msg = msg


def average_rmsd(runs: list[rms.RMSD], n_frames: int) -> list[int]:
    if len(runs) < 1:
        return []

    avg_rmsd = []
    for frame in range(0, n_frames):
        avg_rmsd.append(sum([run.results.rmsd[frame][2] for run in runs]) / len(runs))

    return avg_rmsd


def plot_rmsd(
    ax: Axes,
    rmsd_by_neigh: list[list[int]],
    str_depths: list[str],
    n_frames: int,
):
    if not rmsd_by_neigh or not str_depths:
        return

    for rmsd in rmsd_by_neigh:
        ax.plot(
            [x + 1 for x in range(0, n_frames)],
            rmsd,
        )

    ax.legend(str_depths, title="Neighbourhood depth")
    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD (Å)")


def serialize_rmsd(rmsd: list[int], n_frames: int, filename: str):
    if not rmsd or n_frames < 1:
        return

    with open(filename, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["Frame", "RMSD"])
        for frame in range(0, n_frames):
            writer.writerow([frame + 1, rmsd[frame]])


def compute_rmsf(topology: str, trajectory: str, selection_str: str):
    if topology == "" or trajectory == "":
        return None

    u = mda.Universe(topology, trajectory)

    average = align.AverageStructure(
        u, u, select="protein and name CA", ref_frame=0
    ).run()

    ref = average.results.universe
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()

    c_alphas = u.select_atoms(f"protein and name CA and ({selection_str})")

    rmsf = rms.RMSF(c_alphas).run()

    return rmsf, c_alphas


def serialize_rmsf(rmsf: tuple[list[float], mda.AtomGroup], filename: str):
    with open(filename, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["resi", "RMSF"])
        for resi, rmsf_val in zip(rmsf[1].resids, rmsf[0]):
            writer.writerow([resi, rmsf_val])


def serialize_comparison(
    duration: float,
    full_duration: float,
    n_atoms: int,
    full_n_atoms: int,
    rmsf: tuple[list[float], mda.AtomGroup],
    rmsd: list[int],
    filename: str,
):
    with open(filename, "w") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "Partial sim. duration",
                "Full sim. Duration",
                "Duration change",
                "Partial size",
                "Full size",
                "Size change",
                "Avg RMSD",
                "Avg RMSF",
            ]
        )
        avg_rmsd = sum(rmsd) / len(rmsd)
        avg_rmsf = sum(rmsf[0]) / len(rmsf[0])
        writer.writerow(
            [
                duration,
                full_duration,
                duration / full_duration,
                n_atoms,
                full_n_atoms,
                n_atoms / full_n_atoms,
                avg_rmsd,
                avg_rmsf,
            ]
        )


def average_rmsf(
    runs: list[tuple[rms.RMSF, mda.AtomGroup]],
) -> tuple[list[float], mda.AtomGroup]:
    if not runs:
        return [], mda.AtomGroup()

    rmsf_tot = np.zeros_like(runs[0][0].results.rmsf)
    for rmsf, _ in runs:
        rmsf_tot += np.array(rmsf.results.rmsf)

    return rmsf_tot / len(runs), runs[0][1]


def plot_rmsf(
    ax: Axes,
    rmsf_by_neigh: list[tuple[list[float], mda.AtomGroup]],
    str_neighs: list[str],
):
    if not rmsf_by_neigh or not str_neighs:
        return

    x = np.arange(len(rmsf_by_neigh[0][1].resids))
    bar_width = 0.05
    multiplier = 2
    for (
        rmsf,
        _,
    ), neigh_label in zip(rmsf_by_neigh, str_neighs):
        offset = bar_width * multiplier
        rects = ax.bar(x + offset, rmsf, width=bar_width, label=neigh_label)
        ax.bar_label(rects, fmt="%.2f")
        multiplier += 1

    ax.set_xticks(x + bar_width, rmsf_by_neigh[0][1].resids)
    ax.legend(
        loc="upper left",
        ncols=len(str_neighs),
    )

    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")


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


def plot_duration(ax: Axes, entries: list[int], neighs: list[str]):
    if not entries or not neighs:
        return

    ax.bar(
        neighs,
        entries,
        label=neighs,
    )
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Neighbourhood depth")
    ax.set_ylabel("Duration (s)")


def process_sim(
    molecule_name,
    full_sim,
    minimized,
    sim_dir: str,
    radius: int,
    depth: int,
    paratope_selection: str,
    paratope_hc_selection: str,
    paratope_lc_selection: str,
    frames: int,
):
    # TODO this code does no averaging anymore
    durations = []
    rmsd_runs = []
    rmsd_static_runs = []
    rmsf_h_runs = []
    rmsf_l_runs = []

    print(f"Processing simulation {sim_dir}")

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

    rmsd_runs.append(rms.RMSD(partial, full_sim, select=paratope_selection).run())
    rmsd_static_runs.append(
        rms.RMSD(partial, minimized, select=paratope_selection).run(),
    )

    rmsf_h_runs.append(
        compute_rmsf(
            final_molecule_file,
            os.path.join(sim_dir, "trajectory.dcd"),
            paratope_hc_selection,
        ),  # heavy chain RMSF
    )
    rmsf_l_runs.append(
        compute_rmsf(
            final_molecule_file,
            os.path.join(sim_dir, "trajectory.dcd"),
            paratope_lc_selection,
        ),  # heavy chain RMSF
    )

    # compute averages over multiple runs
    avg_duration = sum(durations) / max(len(durations), 1)
    avg_rmsd = average_rmsd(rmsd_runs, frames)
    avg_rmsd_static = average_rmsd(rmsd_static_runs, frames)
    avg_rmsf_h = average_rmsf(rmsf_h_runs)
    avg_rmsf_l = average_rmsf(rmsf_l_runs)

    # serialize the results
    serialize_rmsd(avg_rmsd, frames, os.path.join(sim_dir, "rmsd.csv"))
    serialize_rmsd(avg_rmsd_static, frames, os.path.join(sim_dir, "rmsd_static.csv"))
    serialize_rmsf(avg_rmsf_h, os.path.join(sim_dir, "rmsf_h.csv"))
    serialize_rmsf(avg_rmsf_l, os.path.join(sim_dir, "rmsf_l.csv"))

    sim_name = f"sim_r{radius}d{depth}"
    cmd.load(
        os.path.join(sim_dir, f"{molecule_name}_minimized_sliced_fixed.pdb"), sim_name
    )
    cmd.load(os.path.join(sim_dir, "trajectory.dcd"), sim_name)

    return (
        avg_duration,
        avg_rmsd,
        avg_rmsd_static,
        avg_rmsf_h,
        avg_rmsf_l,
        partial.atoms.n_atoms,
    )


def compare_sims(base_dir: str, molecule_name: str, n_frames_str: str):
    n_frames = int(n_frames_str)
    molecule_dir = os.path.join(base_dir, molecule_name)

    avg_durations = []
    avg_rmsd_by_neigh = []
    avg_rmsd_by_neigh_static = []
    avg_rmsf_h_by_neigh = []
    avg_rmsf_l_by_neigh = []
    str_neighs = ["Whole molecule"]

    # We have one directory for each simulation
    dirs = [
        f.path for f in os.scandir(molecule_dir) if f.is_dir() and f.name != "ignore"
    ]
    if not dirs:
        print("No simulation directories found")
        return

    # Find the directory containing the full simulation
    full_sim_regex = re.compile(r".*_full")
    full_sim_dir = None
    for dir in dirs:
        if full_sim_regex.match(dir):
            full_sim_dir = dir
            break

    if not full_sim_dir:
        print("No full simulation directory found")
        return

    # Load the full simulation
    final_full_molecule_file = os.path.join(
        molecule_dir, full_sim_dir, f"{molecule_name}_minimized.pdb"
    )
    final_full_traj_file = os.path.join(molecule_dir, full_sim_dir, "trajectory.dcd")
    full_sim = mda.Universe(
        final_full_molecule_file,
        final_full_traj_file,
    )
    minimized = mda.Universe(final_full_molecule_file)

    partial_sim_regex = re.compile(r".*_r([0-9]+)d([0-9]+)")
    tmp = None
    for dir in dirs:
        if partial_sim_regex.match(dir):
            tmp = dir
            break

    if not tmp:
        print("No partial simulation directory found")
        return

    # Load the paratope residues from the json file
    with open(os.path.join(molecule_dir, tmp, "simulated_residues.json")) as f:
        sim_metadata = json.load(f)
        paratope_residues = sim_metadata["paratope"]

    paratope_hc_selection = (
        "name CA and ("
        + " or ".join(
            [
                f"(resid {res} and chainID {chain})"
                for res, chain in paratope_residues
                if chain == "H"
            ]
        )
        + ")"
    )
    paratope_lc_selection = (
        "name CA and ("
        + " or ".join(
            [
                f"(resid {res} and chainID {chain})"
                for res, chain in paratope_residues
                if chain == "L"
            ]
        )
        + ")"
    )

    paratope_selection = (
        "name CA and ("
        + " or ".join(
            [f"(resid {res} and chainID {chain})" for res, chain in paratope_residues]
        )
        + ")"
    )

    # compute the rmsf for the paratope on the full simulation
    rmsf_h_full = compute_rmsf(
        final_full_molecule_file,
        final_full_traj_file,
        paratope_hc_selection,
    )
    if not rmsf_h_full:
        return
    avg_rmsf_h_by_neigh.append(average_rmsf([rmsf_h_full]))

    rmsf_l_full = compute_rmsf(
        final_full_molecule_file,
        final_full_traj_file,
        paratope_lc_selection,
    )
    if not rmsf_l_full:
        return
    avg_rmsf_l_by_neigh.append(average_rmsf([rmsf_l_full]))

    # compute the rmsd for the paratope on the full simulation
    rmsd_static_full = rms.RMSD(full_sim, minimized, select=paratope_selection).run()
    avg_rmsd_by_neigh_static.append(
        list(map(lambda x: x[2], rmsd_static_full.results.rmsd))
    )

    # get the time needed for the full simulation
    with open(os.path.join(molecule_dir, full_sim_dir, "sim_state.csv")) as f:
        final_line = f.readlines()[-1].strip()
        avg_durations.append(float(final_line.split(",")[3]))

    cmd.load(final_full_molecule_file, "full_sim")
    cmd.load(final_full_traj_file, "full_sim")

    for i, dir in enumerate(dirs):
        if partial_sim_regex.match(dir):
            radius = partial_sim_regex.match(dir).group(1)
            depth = partial_sim_regex.match(dir).group(2)
        else:
            print(f"Skipping {dir}")
            continue

        try:
            (
                avg_duration,
                avg_rmsd,
                avg_rmsd_static,
                avg_rmsf_h,
                avg_rmsf_l,
                n_atoms,
            ) = process_sim(
                molecule_name,
                full_sim,
                minimized,
                dir,
                int(radius),
                int(depth),
                paratope_selection,
                paratope_hc_selection,
                paratope_lc_selection,
                n_frames,
            )
        except NoSimulationsException as e:
            print(e.msg)
            continue

        str_neighs.append(f"r{radius}d{depth}")

        avg_durations.append(avg_duration)
        avg_rmsd_by_neigh.append(avg_rmsd)
        avg_rmsd_by_neigh_static.append(avg_rmsd_static)
        avg_rmsf_h_by_neigh.append(avg_rmsf_h)
        avg_rmsf_l_by_neigh.append(avg_rmsf_l)

        serialize_comparison(
            avg_duration,
            avg_durations[0],
            n_atoms,
            full_sim.atoms.n_atoms,
            avg_rmsf_h,
            avg_rmsd,
            os.path.join(dir, "comparison.csv"),
        )

    _, (rmsd_ax, rmsd_static_ax) = plt.subplots(2, 1)
    rmsd_ax.set_title("Partial VS full simulation")
    rmsd_static_ax.set_title("Simulation VS static molecule")
    # compare each frame with the static paratope
    plot_rmsd(rmsd_static_ax, avg_rmsd_by_neigh_static, str_neighs, n_frames)
    # compare each frame of partial simulations with those of the full simulation
    plot_rmsd(rmsd_ax, avg_rmsd_by_neigh, str_neighs[1:], n_frames)

    _, (dur_ax, dur_with_full_ax) = plt.subplots(2, 1)
    dur_ax.set_title("Simulation duration")
    # plot durations only for partial simulations since they are much shorter than the full one
    plot_duration(dur_ax, avg_durations[1:], str_neighs[1:])
    # include the full simulation in the runtime plot
    plot_duration(dur_with_full_ax, avg_durations, str_neighs)

    fig_rmsf_h, rmsf_ax_h = plt.subplots(layout="constrained")
    rmsf_ax_h.set_title("RMSF of the paratope on H")
    # ax_rmsf_h = fig_rmsf_h.add_axes((0.7, 0.05, 0.1, 0.075))
    # btn_rmsf_h = Button(ax_rmsf_h, "Annotate")
    # btn_rmsf_h.on_clicked(lambda _: toggle_annotation(rmsf_ax_h, avg_rmsf_h_by_neigh))
    # buttons.append(btn_rmsf_h)
    plot_rmsf(rmsf_ax_h, avg_rmsf_h_by_neigh, str_neighs)

    fig_rmsf_l, rmsf_ax_l = plt.subplots(layout="constrained")
    rmsf_ax_l.set_title("RMSF of the paratope on L")
    # ax_rmsf_l = fig_rmsf_l.add_axes((0.7, 0.05, 0.1, 0.075))
    # btn_rmsf_l = Button(ax_rmsf_l, "Annotate")
    # btn_rmsf_l.on_clicked(lambda _: toggle_annotation(rmsf_ax_l, avg_rmsf_l_by_neigh))
    # buttons.append(btn_rmsf_l)
    plot_rmsf(rmsf_ax_l, avg_rmsf_l_by_neigh, str_neighs)

    # _, rmsf_ax_full = plt.subplots()
    # rmsf_ax_full.set_title("RMSF on the full molecule")
    # plot_rmsf(rmsf_ax_full, rmsf_full, ["Full molecule"])

    plt.show()


cmd.extend("compare_sims", compare_sims)
