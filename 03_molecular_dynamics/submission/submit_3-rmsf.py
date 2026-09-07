#!/usr/bin/env python3
import os
import glob
import argparse
import shutil

parser = argparse.ArgumentParser(
    description="RMSF submission Script (SLURM-aware)"
)

# -------------------------------------------------------------------------
# REQUIRED INPUTS
# -------------------------------------------------------------------------
parser.add_argument("--top", type=str, default = 'system/popc_0.1mCaCl2_seqfixed_180Apadding/step5_input.psf',
                    help="Topology (psf/pdb)")

parser.add_argument("--rep", type=int, default=0,
                    help="Replica number")

parser.add_argument("--sysname", type=str, default = 'popc_0.1mCaCl2_seqfixed_180Apadding',
                    help="popc_0.1mCaCl2_seqfixed_180Apadding")

parser.add_argument('--trajs', default = None, nargs = '+', help = 'List of TRAJs')

parser.add_argument('--type_trajs', default = '.dcd', help = 'type of trajs')

parser.add_argument('--inpdir', default = None, help = 'Input directory of TRAJs')

parser.add_argument("--stride", type=int, default=1,
                    help="Frame stride")

parser.add_argument("--prefix", default=None,
                    help="Output prefix")

parser.add_argument("--job_name", type=str, default=None)

args = parser.parse_args()


# -------------------------------------------------------------------------
# PATH SETUP
# -------------------------------------------------------------------------
root = args.inpdir
cwd = os.getcwd()
sc = os.path.abspath(os.path.join(cwd, ".."))
os.makedirs(f"{sc}/jobs", exist_ok=True)

job_name = args.job_name if args.job_name else f"RMSF_Rep-{args.rep}"
log_fn = os.path.join(sc, "jobs", f"{job_name}.log")


# -------------------------------------------------------------------------
# FIND TRAJECTORIES
# -------------------------------------------------------------------------
if not args.trajs:
    if args.type_trajs == '.dcd':
        traj_pattern = os.path.join(
            root,
            "simulations",
            args.sysname,
            f"Rep-{args.rep}",
            "*/*.dcd"
            )
    else:
        traj_pattern = os.path.join(
            args.inpdir,
            args.sysname,
            f"Rep-{args.rep}",
            "whole/*.xtc"
            )
    trajs = sorted(glob.glob(traj_pattern))
else:
    trajs = args.trajs

if len(trajs) == 0:
    raise RuntimeError(f"No trajectories found: {traj_pattern}")

print(f"[INFO] Found {len(trajs)} trajectories")


# -------------------------------------------------------------------------
# BUILD RMSF COMMAND
# -------------------------------------------------------------------------
rmsf_cmd = (
    f"python {sc}/scripts/_RMSF.py "
    f"--top {args.top} "
    f"--rep {args.rep} "
    f"--sysname {args.sysname} "
    f"--stride {args.stride} "
    f"--prefix {args.prefix} "
    f"--trajs {' '.join(trajs)}"
)


# -------------------------------------------------------------------------
# DETECT SLURM
# -------------------------------------------------------------------------
has_slurm = shutil.which("sbatch") is not None


# -------------------------------------------------------------------------
# GENERATE JOB SCRIPT
# -------------------------------------------------------------------------
if has_slurm:

    print("[INFO] SLURM detected → submitting job")

    script_text = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH -J {job_name}
#SBATCH -o {log_fn}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate openff

cd {root}

echo "===== Job started $(date) ====="

{rmsf_cmd}

echo "===== Job finished $(date) ====="
"""

    script_path = f"{root}/submission/{job_name}.job"

    with open(script_path, "w") as f:
        f.write(script_text)

    os.system(f"sbatch {script_path}")
    print(f"[INFO] Submitted → {script_path}")


# -------------------------------------------------------------------------
# LOCAL (tmux)
# -------------------------------------------------------------------------
else:

    print("[INFO] SLURM not found → running in tmux")

    script_text = f"""#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate openff

cd {sc}

echo "===== Job started $(date) ====="

{rmsf_cmd}

echo "===== Job finished $(date) ====="
"""

    script_path = f"{sc}/submission/{job_name}.sh"

    with open(script_path, "w") as f:
        f.write(script_text)

    os.chmod(script_path, 0o755)

    os.system(
        f'tmux new-session -d -s {job_name} "bash {script_path} 2>&1 | tee {log_fn}"'
    )

    print(f"[INFO] tmux session started → {job_name}")
