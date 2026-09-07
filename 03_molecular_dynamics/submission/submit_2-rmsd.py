#!/usr/bin/env python3
import os
import glob
import argparse
import shutil   # detect SLURM

parser = argparse.ArgumentParser(
    description="Unified RMSD Submission Script (SLURM-aware)"
)

# -------------------------------------------------------------------------
# GENERAL ARGUMENTS
# -------------------------------------------------------------------------
parser.add_argument("--mode", required=True,
                    choices=["per_residue", "whole_system"],
                    help="Choose RMSD mode")

parser.add_argument("--ref", type=str, default=None,
                    help="Reference PDB (optional). If omitted → use first frame")

parser.add_argument("--top", type=str, default="system/popc_0.1mCaCl2_seqfixed_180Apadding/step5_input.psf",
                    help="Topology file (PSF or PDB)")

parser.add_argument("--indir", nargs='+', required=True,
                    help="Input directories (e.g. minimization simulations)")

parser.add_argument("--type_trajs", default = '.dcd',
                    help="Type of trajs")

parser.add_argument("--rep", type=int, default=0,
                    help="Replica number")

parser.add_argument("--selection", type=str, default=None,
                    help="Selection of group atoms for RMSD/RMSF cals")

parser.add_argument("--out_dir", type=str, default="analysis",
                    help="Directory for output")

parser.add_argument("--out_prefix", type=str, default="rmsd",
                    help="Prefix of output filename")

parser.add_argument("--max_workers", type=int, default=None)

parser.add_argument("--job_name", type=str, default=None)

args = parser.parse_args()


# -------------------------------------------------------------------------
# SETUP
# -------------------------------------------------------------------------
cwd = os.getcwd()
os.chdir(os.path.abspath(os.path.join(cwd, "..")))
cwd = os.getcwd()

os.makedirs(f"{cwd}/submission", exist_ok=True)
os.makedirs(f"{cwd}/jobs", exist_ok=True)

job_name = args.job_name if args.job_name else f"{args.out_prefix}_{args.mode}_Rep-{args.rep}"
log_fn = os.path.join(cwd, "jobs", f"{job_name}.log")

max_workers = args.max_workers if args.max_workers else max(1, os.cpu_count() - 1)


# -------------------------------------------------------------------------
# DETERMINE SYSTEM NAME AND TRAJECTORIES
# -------------------------------------------------------------------------
sysname = os.path.basename(os.path.dirname(args.top))

traj_pattern = []

for d in args.indir:
    # Example: minimization/<sysname>/Rep-0/*/*.dcd
    if args.type_trajs == '.dcd':
        p1 = os.path.join(d, sysname, f"Rep-{args.rep}", "*/*.dcd")
    else:
        type_ = os.path.basename(args.top).split('.')[0]
        p1 = os.path.join(d, sysname, f"Rep-{args.rep}", f"{type_}/*.xtc")

    print (p1)
    traj_pattern.append(p1)

# Expand patterns → lists of lists
trajs = []
for p in traj_pattern:
    trajs.extend(glob.glob(p))

# Sort them
trajs = sorted(trajs)

if len(trajs) == 0:
    raise RuntimeError(f"No trajectories found under: {args.indir}")

print(f"[INFO] Found {len(trajs)} trajectories")

"""
selection = ('protein and backbone and not segid PROH PROR PRON PROP PROJ PROL and '
             'not resid 150:180 and not resid 240:251 and not resid 1:10 and '
             'not (resid 103 and segid PROG PROI PROK PROM PROo PROQ and name O)')
"""

if args.selection is None:
    selection = ('protein and backbone and segid PROA PROB PROC PROD PROE PROF and ' #Residues in the pore
                 'not resid 150:180 and not resid 240:251 and not resid 1:35 and ' 
                 'not resid 208 223')
else:
    selection = args.selection

# -------------------------------------------------------------------------
# OUTPUT FILE
# -------------------------------------------------------------------------
os.makedirs(f"{args.out_dir}/{sysname}/structural/RMSD", exist_ok=True)
out_fn = os.path.join(args.out_dir, sysname, 'structural/RMSD', f"{args.out_prefix}_Rep-{args.rep}")

if args.mode == "per_residue":
    out_fn += ".h5"
else:
    out_fn += ".npy"

# -------------------------------------------------------------------------
# BUILD RMSD COMMAND
# -------------------------------------------------------------------------
rmsd_cmd = f"python {cwd}/scripts/_RMSD.py {args.mode} "

# Only include --ref if user gave one
if args.ref not in [None, "", "None", "none", "NONE"]:
    rmsd_cmd += f"--ref {args.ref} "

rmsd_cmd += (
    f"--top {args.top} "
    f"--traj {' '.join(trajs)} "
    f"--out {out_fn} "
    f"--selection \"{selection}\" "
    )

# Add workers only for per_residue
if args.mode == "per_residue":
    rmsd_cmd += f"--max_workers {max_workers} "


# -------------------------------------------------------------------------
# DETECT SLURM
# -------------------------------------------------------------------------
has_slurm = shutil.which("sbatch") is not None


# -------------------------------------------------------------------------
# GENERATE SCRIPT
# -------------------------------------------------------------------------
if has_slurm:
    print("[INFO] SLURM detected → generating SLURM job")

    wall_time = "10:00:00"
    script_text = f"""#!/bin/bash
#SBATCH -A bio250388
#SBATCH --nodes=1
#SBATCH --ntasks={max_workers}
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time={wall_time}
#SBATCH -J {job_name}
#SBATCH -o {log_fn}
#SBATCH -p shared

module load anaconda/2021.05-py38
conda activate openff

echo "===== Job started at $(date '+%Y-%m-%d, %H:%M:%S') ====="
echo "Using {max_workers} workers"

cd {cwd}

{rmsd_cmd}

echo "===== Job finished at $(date '+%Y-%m-%d, %H:%M:%S') ====="
"""

    script_path = f"{cwd}/submission/{job_name}.job"
    with open(script_path, "w") as f:
        f.write(script_text)

    os.system(f"sbatch {script_path}")
    print(f"[INFO] SLURM job submitted → {script_path}")

else:
    print("[INFO] SLURM not found → running via tmux")

    script_path = f"{cwd}/submission/{job_name}.sh"
    script_text = f"""#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate openff

echo "===== Job started at $(date '+%Y-%m-%d, %H:%M:%S') ====="

cd {cwd}
{rmsd_cmd}

echo "===== Job finished at $(date '+%Y-%m-%d, %H:%M:%S') ====="
"""

    with open(script_path, "w") as f:
        f.write(script_text)

    os.chmod(script_path, 0o755)

    tmux_cmd = f'tmux new-session -d -s {job_name} "bash {script_path} 2>&1 | tee {log_fn}"'
    os.system(tmux_cmd)

    print(f"[INFO] tmux session started → {job_name}")
    print(f"[INFO] Log file: {log_fn}")

