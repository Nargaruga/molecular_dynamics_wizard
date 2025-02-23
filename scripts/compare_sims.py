import os
import os.path
import json

import numpy as np
from matplotlib.axes import Axes
import matplotlib.pyplot as plt

from pymol import cmd
import MDAnalysis as mda  # TODO cite the use of rms and align as required by (https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html)
from MDAnalysis.analysis import rms, align


def compute_rmsd_by_frame(
    mobile: str, target: str, frames: int, is_static: bool = False
):
    rmsd = []
    for frame in range(1, frames + 1):
        rmsd.append(
            cmd.align(
                mobile,
                target,
                cycles=0,
                mobile_state=frame,
                target_state=1 if is_static else frame,
            )[0]
        )
    return rmsd


def average_rmsd_by_frame(runs: list[list[int]], n_frames: int) -> list[int]:
    if len(runs) < 1:
        return []

    avg_rmsd = []
    for frame in range(0, n_frames):
        avg_rmsd.append(sum([run[frame] for run in runs]) / len(runs))

    return avg_rmsd


def plot_rmsd(
    ax: Axes,
    rmsd_by_frame_by_depth: list[list[int]],
    str_depths: list[str],
    n_frames: int,
):
    for rmsd_by_frame in rmsd_by_frame_by_depth:
        ax.plot(
            [x + 1 for x in range(0, n_frames)],
            rmsd_by_frame,
        )
    ax.legend(str_depths, title="Neighbourhood depth")
    ax.set_xlabel("Frame")
    ax.set_ylabel("RMSD (Å)")


def compute_rmsf(topology, trajectory, residues):
    u = mda.Universe(topology, trajectory)

    average = align.AverageStructure(
        u, u, select="protein and name CA", ref_frame=0
    ).run()

    ref = average.results.universe
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()

    selection_str = " or ".join(
        [f"(resid {res} and chainID {chain})" for res, chain in residues]
    )
    c_alphas = u.select_atoms(f"protein and name CA and ({selection_str})")

    rmsf = rms.RMSF(c_alphas).run()

    return rmsf, c_alphas


def average_rmsf(
    runs: list[tuple[rms.RMSF, mda.AtomGroup]],
) -> tuple[list[float], mda.AtomGroup]:
    rmsf_tot = np.zeros_like(runs[0][0].results.rmsf)
    for rmsf, _ in runs:
        rmsf_tot += np.array(rmsf.results.rmsf)

    return rmsf_tot / len(runs), runs[0][1]


def plot_rmsf(
    ax: Axes,
    rmsf_by_depth: list[tuple[list[float], mda.AtomGroup]],
    str_depths: list[str],
):
    for rmsf, c_alphas in rmsf_by_depth:
        ax.plot(c_alphas.resids, rmsf, linestyle="--", marker="o")
        for xy in zip(c_alphas.resids, rmsf):
            ax.annotate("(%s, %.2f)" % xy, xy=xy, textcoords="data")

    ax.legend(str_depths, title="Neighbourhood depth")
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")


def plot_duration(ax: Axes, entries: list[int], depths: list[int]):
    ax.bar(
        depths,
        entries,
        label=depths,
    )
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Neighbourhood depth")
    ax.set_ylabel("Duration (s)")


def process_directory(
    sel_on_minimized: str,
    dir: str,
    depth: int,
    paratope_residues: list[tuple[int, str]],
    frames: int,
):
    print(f"Processing simulations of depth {depth}")

    durations = []
    rmsd_runs = []
    rmsd_static_runs = []
    rmsf_h_runs = []
    rmsf_l_runs = []

    partial_sim = f"partial_sim_d{depth}"
    sel_on_partial = f"paratope_d{depth}"

    paratope_hc_residues = [
        (resi, chain) for (resi, chain) in paratope_residues if chain == "H"
    ]
    paratope_lc_residues = [
        (resi, chain) for (resi, chain) in paratope_residues if chain == "L"
    ]

    # process all simulations performed with this neighbourhood depth
    sim_dirs = [f.path for f in os.scandir(dir) if f.is_dir() and f.name != "ignore"]
    for sim_dir in sim_dirs:
        print(f"Processing simulation {sim_dir}")

        with open(os.path.join(sim_dir, "sim_state.csv")) as f:
            final_line = f.readlines()[-1].strip()
            durations.append(float(final_line.split(",")[3]))

        # load partial sim trajectory
        cmd.load(os.path.join(sim_dir, "ready.pdb"), partial_sim)
        cmd.load_traj(os.path.join(sim_dir, "trajectory.dcd"), partial_sim)

        # select the paratope
        for residue in paratope_residues:
            cmd.select(
                name=sel_on_partial,
                selection=f"byres {partial_sim} and resi {residue[0]} and chain {residue[1]}",
                merge=1,
            )

        rmsd_runs.append(
            compute_rmsd_by_frame(sel_on_partial, sel_on_minimized, frames)
        )
        rmsd_static_runs.append(
            compute_rmsd_by_frame(sel_on_partial, sel_on_minimized, frames, True)
        )

        rmsf_h_runs.append(
            compute_rmsf(
                os.path.join(sim_dir, "ready.pdb"),
                os.path.join(sim_dir, "trajectory.dcd"),
                paratope_hc_residues,
            ),  # heavy chain RMSF
        )
        rmsf_l_runs.append(
            compute_rmsf(
                os.path.join(sim_dir, "ready.pdb"),
                os.path.join(sim_dir, "trajectory.dcd"),
                paratope_lc_residues,
            ),  # heavy chain RMSF
        )

        # cmd.delete(sel_on_partial)

    # cmd.delete(partial_sim)

    # compute averages over multiple runs
    avg_duration = sum(durations) / len(durations)
    avg_rmsd = average_rmsd_by_frame(rmsd_runs, frames)
    avg_rmsd_static = average_rmsd_by_frame(rmsd_static_runs, frames)
    avg_rmsf_h = average_rmsf(rmsf_h_runs)
    avg_rmsf_l = average_rmsf(rmsf_l_runs)

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
    avg_rmsd_by_depth_by_frame = []
    avg_rmsd_by_depth_by_frame_static = []
    avg_rmsf_h_by_depth = []
    avg_rmsf_l_by_depth = []
    # rmsf_full = []
    str_depths = []

    full_sim = "full_sim"
    cmd.load(os.path.join(molecule_dir, "minimized.pdb"), full_sim)
    cmd.load_traj(os.path.join(molecule_dir, "full_trajectory.dcd"), full_sim)

    minimized = "minimized"
    cmd.load(os.path.join(molecule_dir, "minimized.pdb"), minimized)

    with open(os.path.join(molecule_dir, "paratope.json")) as f:
        sim_metadata = json.load(f)
        paratope_residues = sim_metadata["paratope"]

    sel_on_full_sim = "paratope_on_full_sim"
    sel_on_minimized = "paratope_on_minimized"
    for residue in paratope_residues:
        cmd.select(
            name=sel_on_full_sim,
            selection=f"{full_sim} and resi {residue[0]} and chain {residue[1]}",
            merge=1,
        )

        cmd.select(
            name=sel_on_minimized,
            selection=f"{minimized} and resi {residue[0]} and chain {residue[1]}",
            merge=1,
        )

    # compute the rmsd and simulation duration for the various partial simulations
    dirs = [
        f.path for f in os.scandir(molecule_dir) if f.is_dir() and f.name != "ignore"
    ]
    # we assume that the directory name is the neighbourhood depth used for simulations
    str_depths += [os.path.basename(dir) for dir in dirs]
    for i, dir in enumerate(dirs):
        (
            avg_duration,
            avg_rmsd_by_frame,
            avg_rmsd_by_frame_static,
            avg_rmsf_h,
            avg_rmsf_l,
        ) = process_directory(
            sel_on_full_sim,
            dir,
            int(str_depths[i]),
            paratope_residues,
            n_frames,
        )
        avg_durations.append(avg_duration)
        avg_rmsd_by_depth_by_frame.append(avg_rmsd_by_frame)
        avg_rmsd_by_depth_by_frame_static.append(avg_rmsd_by_frame_static)
        avg_rmsf_h_by_depth.append(avg_rmsf_h)
        avg_rmsf_l_by_depth.append(avg_rmsf_l)

    # compute the rmsd for the paratope on the full simulation
    # we need a copy because we can't align an object to itself
    rmsd_by_frame_static_full = compute_rmsd_by_frame(
        sel_on_full_sim, sel_on_minimized, n_frames, True
    )
    avg_rmsd_by_depth_by_frame_static.append(rmsd_by_frame_static_full)

    # compute the rmsf for the paratope on the full simulation
    # rmsf_full.append(compute_rmsf(
    #     os.path.join(molecule_dir, "minimized.pdb"),
    #     os.path.join(molecule_dir, "full_trajectory.dcd"),
    #     paratope_residues,
    # ))

    # get the time needed for the full simulation
    with open(os.path.join(molecule_dir, "full_sim_state.csv")) as f:
        final_line = f.readlines()[-1].strip()
        avg_durations.append(float(final_line.split(",")[3]))

    str_depths.append("Whole molecule")

    _, (rmsd_ax, rmsd_static_ax) = plt.subplots(2, 1)
    rmsd_ax.set_title("Partial VS full simulation")
    rmsd_static_ax.set_title("Simulation VS static molecule")
    # compare each frame with the static paratope
    plot_rmsd(rmsd_static_ax, avg_rmsd_by_depth_by_frame_static, str_depths, n_frames)
    # compare each frame of partial simulations with those of the full simulation
    plot_rmsd(rmsd_ax, avg_rmsd_by_depth_by_frame, str_depths, n_frames)

    _, (dur_ax, dur_with_full_ax) = plt.subplots(2, 1)
    dur_ax.set_title("Simulation duration")
    # plot durations only for partial simulations since they are much shorter than the full one
    plot_duration(dur_ax, avg_durations[:-1], str_depths[:-1])
    # include the full simulation in the runtime plot
    plot_duration(dur_with_full_ax, avg_durations, str_depths)

    _, rmsf_ax_h = plt.subplots()
    rmsf_ax_h.set_title("RMSF of the paratope on H")
    plot_rmsf(rmsf_ax_h, avg_rmsf_h_by_depth, str_depths)

    _, rmsf_ax_l = plt.subplots()
    rmsf_ax_l.set_title("RMSF of the paratope on L")
    plot_rmsf(rmsf_ax_l, avg_rmsf_l_by_depth, str_depths)

    # _, rmsf_ax_full = plt.subplots()
    # rmsf_ax_full.set_title("RMSF on the full molecule")
    # plot_rmsf(rmsf_ax_full, rmsf_full, ["Full molecule"])

    plt.show()


cmd.extend("compare_sims", compare_sims)
