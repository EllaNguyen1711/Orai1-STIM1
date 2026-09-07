import os
import glob
import time
import argparse

import MDAnalysis as mda
from MDAnalysis.analysis import align

## -- FUNCTIONS -- ##

def align_traj_stride_write(topology, traj, ref, out_traj, out_pdb,
                            selection="protein and name CA",
                            stride=1):
    """
    Align trajectory to a reference with stride, write aligned trajectory
    and a PDB of the aligned first frame.
    """
    u     = mda.Universe(topology, traj)
    ref_u = mda.Universe(topology, ref)

    align.AlignTraj(
        u,
        ref_u,
        select=selection,
        filename=out_traj,
        step=stride,
        in_memory=False
    ).run()

    u.trajectory[0]
    with mda.Writer(out_pdb) as W:
        W.write(u.atoms)


def write_chunked_xtc(topology, trajs, selection="protein",
                      outdir=".", basename="TRAJ",
                      stride=1, chunk_size=2000,
                      continuous=False, write_pdb=True, overwrite=False):
    """
    Write a selection of atoms from one or more input trajectories into a list
    of XTC files, each containing at most `chunk_size` frames.

    Parameters
    ----------
    topology : str
        Topology file (e.g. step5_input.psf).
    trajs : str or list of str
        Input trajectory file(s), e.g. the sorted list of TRAJ*/sim.dcd.
    selection : str
        MDAnalysis atom selection to write (e.g. "protein", "all").
    outdir : str
        Output directory. Created if missing.
    basename : str
        Prefix for the output files.
    stride : int
        Keep every `stride`-th frame. Chunking is applied AFTER striding, so
        each output XTC holds `chunk_size` written frames.
    chunk_size : int
        Number of frames per output XTC (2000 by default).
    continuous : bool
        If False (default), each input trajectory is chunked independently:
            TRAJ000.part0000.xtc, TRAJ000.part0001.xtc, TRAJ001.part0000.xtc, ...
        If True, all input trajectories are concatenated first and chunked as a
        single continuous stream, so every file except the last has exactly
        `chunk_size` frames:
            TRAJ.part0000.xtc, TRAJ.part0001.xtc, ...
    write_pdb : bool
        Also write a PDB of the selection (first frame) per group.
    overwrite : bool
        If False, a group whose first chunk already exists is skipped.

    Returns
    -------
    list of str
        Paths of the XTC files written (or already present when skipped).
    """
    if isinstance(trajs, str):
        trajs = [trajs]
    if not trajs:
        raise ValueError("No input trajectories were given.")

    os.makedirs(outdir, exist_ok=True)

    # One "group" == one Universe == one independent chunk counter.
    if continuous:
        groups = [(basename, list(trajs))]
    else:
        groups = [(f"{basename}{i:03d}", [t]) for i, t in enumerate(trajs)]

    all_written = []

    for name, group in groups:
        first_chunk = os.path.join(outdir, f"{name}.part0000.xtc")
        if os.path.isfile(first_chunk) and not overwrite:
            print(f"  [skip] {name}: {os.path.basename(first_chunk)} already exists")
            all_written += sorted(glob.glob(os.path.join(outdir, f"{name}.part*.xtc")))
            continue

        u     = mda.Universe(topology, group)
        atoms = u.select_atoms(selection)
        if atoms.n_atoms == 0:
            raise ValueError(f"Selection '{selection}' matched 0 atoms in {topology}")

        # reference PDB of the selection (first frame)
        if write_pdb:
            u.trajectory[0]
            pdbpath = os.path.join(outdir, f"{name}.pdb")
            if overwrite or not os.path.isfile(pdbpath):
                with mda.Writer(pdbpath, atoms.n_atoms) as P:
                    P.write(atoms)

        writer      = None
        outtraj     = None
        n_chunks    = 0
        n_in_chunk  = 0
        n_total     = 0

        try:
            for ts in u.trajectory[::stride]:
                if n_in_chunk == 0:
                    outtraj = os.path.join(outdir, f"{name}.part{n_chunks:04d}.xtc")
                    writer  = mda.Writer(outtraj, atoms.n_atoms)

                writer.write(atoms)
                n_in_chunk += 1
                n_total    += 1

                if n_in_chunk == chunk_size:
                    writer.close()
                    writer = None
                    all_written.append(outtraj)
                    print(f"  [write] {os.path.basename(outtraj)}  ({n_in_chunk} frames)")
                    n_chunks   += 1
                    n_in_chunk  = 0
        finally:
            if writer is not None:
                writer.close()
                all_written.append(outtraj)
                print(f"  [write] {os.path.basename(outtraj)}  ({n_in_chunk} frames, last chunk)")
                n_chunks += 1

        print(f"  {name}: {n_total} frames written into {n_chunks} XTC file(s) "
              f"[{atoms.n_atoms} atoms, stride={stride}]")

    return all_written

## -- END FUNCTIONS -- ##

# selection keyword -> (output subfolder, MDAnalysis selection string)
SELECTIONS = {
    'all':     ('whole',   'all'),
    'protein': ('protein', 'protein'),
    'protmem': ('protmem', 'protein or resname POPC or resname CAL CLA'),
}

parser = argparse.ArgumentParser(description="POSTPROCESSING-MD")
parser.add_argument("--sysname",    type=str, default='popc_0.1mCaCl2_seqfixed_180Apadding')
parser.add_argument("--rep",        choices=['0', '1', '2', '3'], help='Replica for TRAJs')
parser.add_argument("--selection",  type=str, choices=list(SELECTIONS), required=True,
                    help="Which atoms to write out")
parser.add_argument("--outdir",     type=str, help="Output directory")
parser.add_argument("--inpdir",     type=str, required=True,
                    help="Input directory where the simulation of GIVEN SYSTEM was stored at")
parser.add_argument("--mutation",   action=argparse.BooleanOptionalAction, default=True)
parser.add_argument("--stride",     type=int, default=1, help="Write every N-th frame")
parser.add_argument("--chunk_size", type=int, default=2000,
                    help="Number of frames per output XTC file")
parser.add_argument("--continuous", action="store_true",
                    help="Concatenate all input DCDs and chunk as one continuous stream")
parser.add_argument("--overwrite",  action="store_true",
                    help="Rewrite chunks even if they already exist")
args = parser.parse_args()

start = time.time()

# 1. Extract and define input trajectory and topology files
inpdir = os.path.join(os.path.abspath(args.inpdir), args.sysname, f'Rep-{args.rep}')
trajs  = sorted(glob.glob(os.path.join(inpdir, 'TRAJ*/sim.dcd')))
top    = os.path.join(os.path.dirname(os.path.abspath(args.inpdir)),
                      'system', args.sysname, 'step5_input.psf')

if not trajs:
    raise SystemExit(f'No sim.dcd found under {inpdir}')

print('Current simulations collected:')
for t in trajs:
    print('   ', t)

dirs = {}
dirs['analysis'] = os.path.join(os.path.abspath('analysis'), args.sysname)
dirs['postsim']  = os.path.join(os.path.abspath('postprocessed_simulation'), args.sysname)
dirs['model']    = os.path.abspath('model')

if args.mutation:
    ref = os.path.join(dirs['model'], f'{args.sysname}.pdb')
else:
    ref = os.path.join(dirs['model'], 'WTprot.model.pdb')

# 2. Write the selected atoms as 2000-frame XTC chunks
base_outdir      = args.outdir if os.path.exists(args.outdir) else os.path.join(dirs['postsim'], f'Rep-{args.rep}')
subdir, selstring = SELECTIONS[args.selection]
outdir           = os.path.join(base_outdir, subdir)

print(f'\nWriting "{args.selection}" ({selstring}) -> {outdir}')
print(f'stride = {args.stride}, chunk_size = {args.chunk_size}, '
      f'mode = {"continuous" if args.continuous else "per-trajectory"}\n')

written = write_chunked_xtc(
    topology=top,
    trajs=trajs,
    selection=selstring,
    outdir=outdir,
    basename='TRAJ',
    stride=args.stride,
    chunk_size=args.chunk_size,
    continuous=args.continuous,
    write_pdb=True,
    overwrite=args.overwrite,
)

end = time.time()
print(f'\n{len(written)} XTC file(s) total')
print(f'Job done within {(end - start)/60:.2f} minutes')
