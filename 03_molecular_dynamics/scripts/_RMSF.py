import time
import argparse
import numpy as np
import MDAnalysis as mda
import os
from MDAnalysis.analysis import rms, align
from MDAnalysis.transformations import center_in_box, fit_rot_trans
#from MDAnalysis.transformations.nojump import NoJump


parser = argparse.ArgumentParser(description='RMSF for selected TRAJs')

parser.add_argument('--rep', type = str, help = 'Replica of TRAJs')
parser.add_argument('--sysname', type = str, help = 'Name of system')
parser.add_argument('--top', type=str, required=True, help = 'Topology file' )
parser.add_argument('--trajs', default = None, nargs='+', help = 'List of TRAJs')
parser.add_argument('--stride', type=int, default=1, help = 'Stride used for data')
parser.add_argument('--prefix', help = 'PREFIX for output filename')

args = parser.parse_args()

start_time = time.time()

psf = args.top
stride = args.stride


# ========================
# Load trajectories
# ========================
trajs = args.trajs
u = mda.Universe(psf, trajs)

# ========================
# Atom selections
# ========================
protein = u.select_atoms('protein')
mobile = u.select_atoms('protein and name CA')


# ========================
# Average structure for alignment
# ========================
avg = align.AverageStructure(
    u, u,
    select='protein and name CA',
    ref_frame=0
).run(step=stride)

ref = avg.results.universe
ref_mobile = ref.select_atoms('protein and name CA')


# ========================
# Transformations
# ========================
transforms = [
    fit_rot_trans(mobile, ref_mobile),
]

u.trajectory.add_transformations(*transforms)


# ========================
# RMSF
# ========================
outdir = os.path.join('analysis', args.sysname, 'structural', f'Rep-{args.rep}')
os.makedirs(outdir, exist_ok=True)

name = f'{args.prefix}.RMSF_Rep-{args.rep}' if args.prefix not in [None, 'None'] else f'RMSF_Rep-{args.rep}'
rmsf_path = os.path.join(outdir, f'{name}.npy')

rmsf_obj = rms.RMSF(mobile).run(step=stride)
np.save(rmsf_path, rmsf_obj.results.rmsf)

print("Saved:", rmsf_path)

print(f"Finished in {(time.time()-start_time)/60:.2f} min")
