import os
import os.path
import json

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
    rmsd_by_depth: list[list[int]],
    str_depths: list[str],
    n_frames: int,
):
    if not rmsd_by_depth or not str_depths:
        return

    for rmsd in rmsd_by_depth:
        ax.plot(
            [x + 1 for x in range(0, n_frames)],
            rmsd,
        )
    ax.legend(str_depths, title="Neighbourhood depth")
    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD (Å)")


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
    rmsf_by_depth: list[tuple[list[float], mda.AtomGroup]],
    str_depths: list[str],
):
    if not rmsf_by_depth or not str_depths:
        return

    for rmsf, c_alphas in rmsf_by_depth:
        ax.plot(c_alphas.resids, rmsf, linestyle="--", marker="o")

    annotate_rmsf(ax, rmsf_by_depth)

    ax.legend(str_depths, title="Neighbourhood depth")
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")


def annotate_rmsf(ax: Axes, rmsf_by_depth: list[tuple[list[float], mda.AtomGroup]]):
    if not rmsf_by_depth:
        return

    for rmsf, c_alphas in rmsf_by_depth:
        for xy in zip(c_alphas.resids, rmsf):
            ax.annotate("(%s, %.2f)" % xy, xy=xy, textcoords="data")


def remove_annotations(ax: Axes):
    for annotation in ax.texts:
        annotation.remove()


def toggle_annotation(ax: Axes, rmsf_by_depth):
    if ax.texts:
        remove_annotations(ax)
    else:
        annotate_rmsf(ax, rmsf_by_depth)

    plt.draw()


def plot_duration(ax: Axes, entries: list[int], depths: list[int]):
    if not entries or not depths:
        return

    ax.bar(
        depths,
        entries,
        label=depths,
    )
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Neighbourhood depth")
    ax.set_ylabel("Duration (s)")


def process_directory(
    full_sim,
    minimized,
    dir: str,
    depth: int,
    paratope_selection: str,
    paratope_hc_selection: str,
    paratope_lc_selection: str,
    frames: int,
):
    print(f"Processing simulations of depth {depth}")

    durations = []
    rmsd_runs = []
    rmsd_static_runs = []
    rmsf_h_runs = []
    rmsf_l_runs = []

    # process all simulations performed with this neighbourhood depth
    sim_dirs = [f.path for f in os.scandir(dir) if f.is_dir() and f.name != "ignore"]
    if not sim_dirs:
        raise NoSimulationsException(f"No simulations found in {dir}")

    for sim_dir in sim_dirs:
        print(f"Processing simulation {sim_dir}")

        with open(os.path.join(sim_dir, "sim_state.csv")) as f:
            final_line = f.readlines()[-1].strip()
            durations.append(float(final_line.split(",")[3]))

        partial = mda.Universe(
            os.path.join(sim_dir, "ready.pdb"), os.path.join(sim_dir, "trajectory.dcd")
        )

        rmsd_runs.append(rms.RMSD(partial, full_sim, select=paratope_selection).run())
        rmsd_static_runs.append(
            rms.RMSD(partial, minimized, select=paratope_selection).run()
        )

        rmsf_h_runs.append(
            compute_rmsf(
                os.path.join(sim_dir, "ready.pdb"),
                os.path.join(sim_dir, "trajectory.dcd"),
                paratope_hc_selection,
            ),  # heavy chain RMSF
        )
        rmsf_l_runs.append(
            compute_rmsf(
                os.path.join(sim_dir, "ready.pdb"),
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

    cmd.load(os.path.join(sim_dirs[0], "ready.pdb"), f"sim_d{depth}")
    cmd.load(os.path.join(sim_dirs[0], "trajectory.dcd"), f"sim_d{depth}")

    return (
        avg_duration,
        avg_rmsd,
        avg_rmsd_static,
        avg_rmsf_h,
        avg_rmsf_l,
        # rmsf_runs[0],  # TODO use avg instead
    )


def compare_sims(molecule_name: str, n_frames_str: str):
    n_frames = int(n_frames_str)
    base = "/home/leo/anaconda3/envs/md_test/lib/python3.9/site-packages/pymol/wizard/tmp/sim_comparison"
    molecule_dir = os.path.join(base, molecule_name)

    avg_durations = []
    avg_rmsd_by_depth = []
    avg_rmsd_by_depth_static = []
    avg_rmsf_h_by_depth = []
    avg_rmsf_l_by_depth = []
    str_depths = []

    full_sim = minimized = mda.Universe(
        os.path.join(molecule_dir, "minimized.pdb"),
        os.path.join(molecule_dir, "full_trajectory.dcd"),
    )

    minimized = mda.Universe(os.path.join(molecule_dir, "minimized.pdb"))

    with open(os.path.join(molecule_dir, "interaction_zone.json")) as f:
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

    # compute the rmsd and simulation duration for the various partial simulations
    dirs = [
        f.path for f in os.scandir(molecule_dir) if f.is_dir() and f.name != "ignore"
    ]
    if not dirs:
        print("No simulation directories found")
        return

    for i, dir in enumerate(dirs):
        # we assume that the directory name corresponds to the simulation depth
        depth = os.path.basename(dir)
        try:
            (
                avg_duration,
                avg_rmsd,
                avg_rmsd_static,
                avg_rmsf_h,
                avg_rmsf_l,
            ) = process_directory(
                full_sim,
                minimized,
                dir,
                int(depth),
                paratope_selection,
                paratope_hc_selection,
                paratope_lc_selection,
                n_frames,
            )
        except NoSimulationsException as e:
            print(e.msg)
            continue

        str_depths.append(depth)

        avg_durations.append(avg_duration)
        avg_rmsd_by_depth.append(avg_rmsd)
        avg_rmsd_by_depth_static.append(avg_rmsd_static)
        avg_rmsf_h_by_depth.append(avg_rmsf_h)
        avg_rmsf_l_by_depth.append(avg_rmsf_l)

    rmsd_static_full = rms.RMSD(full_sim, minimized, select=paratope_selection).run()
    avg_rmsd_by_depth_static.append(average_rmsd([rmsd_static_full], n_frames))

    # compute the rmsf for the paratope on the full simulation
    rmsf_h_full = compute_rmsf(
        os.path.join(molecule_dir, "minimized.pdb"),
        os.path.join(molecule_dir, "full_trajectory.dcd"),
        paratope_hc_selection,
    )
    if not rmsf_h_full:
        return
    avg_rmsf_h_by_depth.append(average_rmsf([rmsf_h_full]))

    rmsf_l_full = compute_rmsf(
        os.path.join(molecule_dir, "minimized.pdb"),
        os.path.join(molecule_dir, "full_trajectory.dcd"),
        paratope_lc_selection,
    )
    if not rmsf_l_full:
        return
    avg_rmsf_l_by_depth.append(average_rmsf([rmsf_l_full]))

    # get the time needed for the full simulation
    with open(os.path.join(molecule_dir, "full_sim_state.csv")) as f:
        final_line = f.readlines()[-1].strip()
        avg_durations.append(float(final_line.split(",")[3]))

    str_depths.append("Whole molecule")

    cmd.load(os.path.join(molecule_dir, "minimized.pdb"), "full_sim")
    cmd.load(os.path.join(molecule_dir, "full_trajectory.dcd"), "full_sim")

    _, (rmsd_ax, rmsd_static_ax) = plt.subplots(2, 1)
    rmsd_ax.set_title("Partial VS full simulation")
    rmsd_static_ax.set_title("Simulation VS static molecule")
    # compare each frame with the static paratope
    plot_rmsd(rmsd_static_ax, avg_rmsd_by_depth_static, str_depths, n_frames)
    # compare each frame of partial simulations with those of the full simulation
    plot_rmsd(rmsd_ax, avg_rmsd_by_depth, str_depths, n_frames)

    _, (dur_ax, dur_with_full_ax) = plt.subplots(2, 1)
    dur_ax.set_title("Simulation duration")
    # plot durations only for partial simulations since they are much shorter than the full one
    plot_duration(dur_ax, avg_durations[:-1], str_depths[:-1])
    # include the full simulation in the runtime plot
    plot_duration(dur_with_full_ax, avg_durations, str_depths)

    fig_rmsf_h, rmsf_ax_h = plt.subplots()
    fig_rmsf_h.subplots_adjust(bottom=0.2)
    rmsf_ax_h.set_title("RMSF of the paratope on H")
    ax_rmsf_h = fig_rmsf_h.add_axes((0.7, 0.05, 0.1, 0.075))
    btn_rmsf_h = Button(ax_rmsf_h, "Annotate")
    btn_rmsf_h.on_clicked(lambda _: toggle_annotation(rmsf_ax_h, avg_rmsf_h_by_depth))
    buttons.append(btn_rmsf_h)
    plot_rmsf(rmsf_ax_h, avg_rmsf_h_by_depth, str_depths)

    fig_rmsf_l, rmsf_ax_l = plt.subplots()
    fig_rmsf_l.subplots_adjust(bottom=0.2)
    rmsf_ax_l.set_title("RMSF of the paratope on L")
    ax_rmsf_l = fig_rmsf_l.add_axes((0.7, 0.05, 0.1, 0.075))
    btn_rmsf_l = Button(ax_rmsf_l, "Annotate")
    btn_rmsf_l.on_clicked(lambda _: toggle_annotation(rmsf_ax_l, avg_rmsf_l_by_depth))
    buttons.append(btn_rmsf_l)
    plot_rmsf(rmsf_ax_l, avg_rmsf_l_by_depth, str_depths)

    # _, rmsf_ax_full = plt.subplots()
    # rmsf_ax_full.set_title("RMSF on the full molecule")
    # plot_rmsf(rmsf_ax_full, rmsf_full, ["Full molecule"])

    plt.show()


cmd.extend("compare_sims", compare_sims)
