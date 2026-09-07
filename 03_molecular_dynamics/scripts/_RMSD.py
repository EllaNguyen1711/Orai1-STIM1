#!/usr/bin/env python3
"""
RMSD for whole system and per residue!
Author: Ella
"""

import os
import argparse
import numpy as np
import time
import h5py as h5
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms
#from MDAnalysis.transformations.nojump import NoJump
#from MDAnalysis.transformations import center_in_box, fit_rot_trans
from concurrent.futures import ProcessPoolExecutor

def get_reference_universe(ref, top, traj):
    """
    Return a reference Universe.

    If ref is provided → load it normally.
    If ref is None → use first frame of trajectory as reference.
    """
    if ref is not None:
        print(f"[INFO] Using provided reference structure: {ref}")
        return mda.Universe(ref)

    print("[INFO] No reference provided → using FIRST FRAME of trajectory")
    # Load system including trajectory
    u = mda.Universe(top, traj)
    refu = mda.Universe(top)
    refu.load_new(u.trajectory[0].positions.copy())
    return refu

def compute_whole_system_rmsd(
    ref,
    top,
    trajs,
    selection,
    out_fn,
    align_selection=None
):

    """
    Parameters
    ----------
    ref : str or Universe
        reference pdb or trajectory
    top : str
        topology
    trajs : list[str]
        trajectories
    selection : str
        atoms for RMSD calculation
    out_fn : str
        output .npy file
    align_selection : str
        atoms used for alignment (default = selection)
    """

    print("\n=== WHOLE SYSTEM RMSD (ALIGNED) ===")
    print(f"Topology: {top}")
    print(f"N trajs: {len(trajs)}")
    print(f"RMSD selection: {selection}")

    sim = mda.Universe(top, trajs)

    refu = get_reference_universe(ref, top, trajs[0])

    # default alignment selection
    if align_selection is None:
        align_selection = selection

    print(f"[INFO] Aligning trajectory using: {align_selection}")

    aligner = align.AlignTraj(
        sim,
        refu,
        select=align_selection,
        in_memory=True
    ).run()

    print("[INFO] Computing RMSD...")
    R = rms.RMSD(
        sim,
        refu,
        select=selection
    ).run()

    rmsd_vals = R.results.rmsd[:, -1]

    np.save(out_fn, rmsd_vals)

    print(f"[Saved] {out_fn}")

    return rmsd_vals


def compute_per_residue_rmsd(ref, top, trajs, selection, out_fn, max_workers=None):

    print("\n=== PER-RESIDUE RMSD (NoJump PBC-fixed) ===")
    print(f"Topology: {top}")
    print(f"Trajectories: {len(trajs)}")
    print(f"Selection: {selection}")
    print(f"Output: {out_fn}")

    # ------------------------------------------------------------
    # Load reference and simulation
    # ------------------------------------------------------------
    refu = get_reference_universe(ref, top, trajs[0])
    sim = mda.Universe(top, trajs)

    # ------------------------------------------------------------
    # Apply PBC corrections BEFORE multiprocessing
    # ------------------------------------------------------------
    mobile_all = sim.atoms

    trans = [
        NoJump(mobile_all),           # remove jumps across boundaries
        center_in_box(mobile_all),    # optional but stabilizes RMSD/RMSF
    ]
    sim.trajectory.add_transformations(*trans)

    # ------------------------------------------------------------
    # Selections
    # ------------------------------------------------------------
    ref_sel = refu.select_atoms(selection)
    traj_sel = sim.select_atoms(selection)

    n_res = len(ref_sel.residues)
    n_frames = len(sim.trajectory)

    print("\n[2/5] Precomputing reference coords")
    ref_coords = [res.atoms.positions.copy() for res in ref_sel.residues]

    # ------------------------------------------------------------
    # Extract PBC-corrected coords for each frame BEFORE MP
    # ------------------------------------------------------------
    print("\n[3/5] Extracting PBC-corrected trajectory coords")
    traj_by_res = [[] for _ in range(n_res)]

    for ts in sim.trajectory:
        for j, res in enumerate(traj_sel.residues):
            traj_by_res[j].append(res.atoms.positions.copy())

    traj_by_res = [np.array(frames) for frames in traj_by_res]

    # ------------------------------------------------------------
    # Worker for multiprocessing
    # ------------------------------------------------------------
    def worker(j):
        refc = ref_coords[j]
        trajc = traj_by_res[j]

        if trajc.shape[1:] != refc.shape:
            return np.full(trajc.shape[0], np.nan, dtype=np.float32)

        diff = trajc - refc[None, :, :]
        return np.sqrt((diff**2).sum(axis=2).mean(axis=1)).astype(np.float32)

    maxw = max_workers if max_workers else max(1, os.cpu_count() - 1)

    print(f"\n[4/5] Multiprocessing with {maxw} workers...")
    with ProcessPoolExecutor(max_workers=maxw) as pool:
        results = list(pool.map(worker, range(n_res)))

    rmsd_matrix = np.column_stack(results)

    # ------------------------------------------------------------
    # Save to HDF5
    # ------------------------------------------------------------
    print("\n[5/5] Saving HDF5 results")
    os.makedirs(os.path.dirname(out_fn) or ".", exist_ok=True)

    with h5.File(out_fn, "w") as f:
        f["rmsd_matrix"] = rmsd_matrix
        f["residue_names"] = np.array([r.resname for r in ref_sel.residues], dtype="S10")
        f["residue_ids"] = np.array([r.resid for r in ref_sel.residues])
        f["mean"] = np.nanmean(rmsd_matrix, axis=0)
        f["std"] = np.nanstd(rmsd_matrix, axis=0)

    print("[DONE] Saved per-residue RMSD (NoJump applied)")
    return rmsd_matrix


# ======================================================================
# CLI
# ======================================================================
def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)

    def add_common(p):
        p.add_argument("--ref", required=False, default=None,
                       help="Reference PDB (optional). If omitted → first frame used.")
        p.add_argument("--top", required=True)
        p.add_argument("--traj", nargs="+", required=True)
        p.add_argument("--selection", default="protein and backbone")
        p.add_argument("--out", required=True)

    pr = sub.add_parser("per_residue")
    add_common(pr)
    pr.add_argument("--max_workers", type=int, default=None)

    ws = sub.add_parser("whole_system")
    add_common(ws)

    args = parser.parse_args()

    if args.mode == "per_residue":
        compute_per_residue_rmsd(args.ref, args.top, args.traj, args.selection, args.out, args.max_workers)
    else:
        compute_whole_system_rmsd(args.ref, args.top, args.traj, args.selection, args.out)


if __name__ == "__main__":
    main()

