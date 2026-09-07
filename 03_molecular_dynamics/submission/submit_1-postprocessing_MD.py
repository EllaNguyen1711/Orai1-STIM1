#!/usr/bin/env python3
import os
import glob
import argparse
import subprocess

parser = argparse.ArgumentParser(
    description="Submit postprocessing-MD job via tmux session."
)
parser.add_argument('--sys',
                    default='popc_0.1mCaCl2_seqfixed_180Apadding')
parser.add_argument('--rep',
                    choices=['0', '1', '2', '3'],
                    required=True)
parser.add_argument('--selection',
                    default='all',
                    help='Atom selection (e.g. all, protein, backbone). Default: all')
parser.add_argument('--stride',
                    default=1,
                    help='Stride (N) for writing trajs every N frames. Default: 1')
parser.add_argument('--inpdir',
                    required=True,
                    help='Input directory where simulations are stored')
parser.add_argument('--outdir', 
                    default = False,
                    type = str,
                    help='Output directory. Default: postprocessed_simulation')
parser.add_argument('--mutation',
                    action='store_true',
                    default=True,
                    help='Use mutant reference model (default: True)')
parser.add_argument('--no_mutation',
                    dest='mutation',
                    action='store_false',
                    help='Use WT reference model instead')
parser.add_argument('--conda_env',
                    default='openff',
                    help='Conda environment to activate. Default: openff')
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
root   = os.path.abspath(os.path.join(os.getcwd(), '..'))
inpdir = os.path.abspath(args.inpdir)
outdir = args.outdir if args.outdir else os.path.join(os.path.abspath('postprocessed_simulation'), args.sys, f'Rep-{args.rep}')
script = os.path.join(root, 'scripts', 'postprocessing_MD.py')

submission_dir = os.path.join(root, 'submission')
jobs_dir       = os.path.join(root, 'jobs')
os.makedirs(submission_dir, exist_ok=True)
os.makedirs(jobs_dir,       exist_ok=True)

# ---------------------------------------------------------------------------
# Sanity check — make sure trajectories exist before submitting
# ---------------------------------------------------------------------------
traj_pattern = os.path.join(inpdir, args.sys, f'Rep-{args.rep}', 'TRAJ*/sim.dcd')
trajs = sorted(glob.glob(traj_pattern))
if not trajs:
    raise RuntimeError(f"No trajectories found matching: {traj_pattern}")
print(f"Trajectories found: {len(trajs)}")

# ---------------------------------------------------------------------------
# Job naming
# ---------------------------------------------------------------------------
job_name = f'postMD_{args.sys}_Rep-{args.rep}_{args.selection}'
log_fn   = os.path.join(jobs_dir,       f'{job_name}.log')
bash_fn  = os.path.join(submission_dir, f'{job_name}.sh')

# ---------------------------------------------------------------------------
# Optional args forwarded to postprocessing_MD.py
# ---------------------------------------------------------------------------
optional_args = [f'--selection {args.selection}']
if not args.mutation:
    optional_args.append('--mutation False')
if not args.stride == 1:
    optional_args.append(f'--stride {args.stride}')

optional_str = " \\\n    ".join(optional_args)

# ---------------------------------------------------------------------------
# Write bash script
# ---------------------------------------------------------------------------
bash_script = f"""#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate {args.conda_env}

echo "===== Job started at $(date) ====="
echo "Host     : $(hostname)"
echo "CPU cores: $(nproc)"
echo "System   : {args.sys}"
echo "Replica  : {args.rep}"
echo "Selection: {args.selection}"
echo "Trajs    : {len(trajs)}"

cd {root}

python -u {script} \\
    --sysname {args.sys} \\
    --rep {args.rep} \\
    --inpdir "{inpdir}" \\
    --outdir "{outdir}" \\
    {optional_str}

echo "===== Job finished at $(date) ====="
"""

with open(bash_fn, 'w') as f:
    f.write(bash_script)
os.chmod(bash_fn, 0o755)

# ---------------------------------------------------------------------------
# Launch in tmux
# ---------------------------------------------------------------------------
subprocess.run(
    f'tmux new-session -d -s {job_name} "bash {bash_fn} 2>&1 | tee {log_fn}"',
    shell=True,
    check=True
)

print(f"Job launched in tmux session : {job_name}")
print(f"Bash script                  : {bash_fn}")
print(f"Log file                     : {log_fn}")
