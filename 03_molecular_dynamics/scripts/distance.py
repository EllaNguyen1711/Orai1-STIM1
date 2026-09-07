import os
import h5py as h5
import numpy as np
import MDAnalysis as mda
from multiprocessing import Pool, cpu_count
import argparse
import sys
sys.path.insert(1, '/media/volume/Orai1_rep0/Orai1_analysis/scripts')
import _Distance

def _compute_chunk(args):
    psf, trajs, selection, start, end, pairs = args
    func = getattr(_Distance, f"dis_{selection}")
    dist = func(psf, trajs, start=start, end=end, pairs = pairs)
    return start, dist

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--psf', required=True)
    parser.add_argument('--traj_fns', nargs='+', required=True)
    parser.add_argument('--selection', required=True)
    parser.add_argument('--list_pairs', default=None, help = 'List of pairs given by a .pickle file')
    parser.add_argument('--out_fn', required=True)
    parser.add_argument('--workers', type=int, default=1,
                        help='Number of parallel workers. 1 = no multiprocessing (default)')
    parser.add_argument('--chunk_size', type=int, default=20000,
                        help='Frames per chunk (only used when --workers > 1)')
    args = parser.parse_args()

    u = mda.Universe(args.psf, args.traj_fns)
    n_frames = len(u.trajectory)
    print(f"[info] frames:     {n_frames}")
    print(f"[info] workers:    {args.workers}")

    if args.workers == 1:
        # --- single process: full trajectory in one call ---
        print("[info] mode: single process")
        func = getattr(_Distance, f"dis_{args.selection}")
        dist = func(args.psf, args.traj_fns, start=0, end=n_frames, pairs = args.list_pairs)
        results = [(0, dist)]

    else:
        # --- multiprocessing: split into chunks ---
        chunk_size = args.chunk_size
        print(f"[info] mode: multiprocessing (chunk_size={chunk_size})")
        jobs = [
            (args.psf, args.traj_fns, args.selection, s, min(s + chunk_size, n_frames), args.list_pairs)
            for s in range(0, n_frames, chunk_size)
        ]
        print(f"[info] total chunks: {len(jobs)}")
        with Pool(args.workers) as pool:
            results = pool.map(_compute_chunk, jobs)

    # --- write output ---
    dist_full = np.empty(n_frames, dtype=np.float32)
    for s, arr in results:
        dist_full[s:s + len(arr)] = arr

    with h5.File(args.out_fn, 'a') as f:
        f.create_dataset("distance", data=dist_full, dtype=np.float32)

    print(f"[DONE] saved → {args.out_fn}")

if __name__ == "__main__":
    main()
