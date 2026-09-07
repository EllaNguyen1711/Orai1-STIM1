import os
import glob
import argparse
import subprocess


parser = argparse.ArgumentParser(
    description="Prepare parallel distance analysis job (chunked multiprocessing)."
)

parser.add_argument('--sys', default='WT_boltz2', help = 'System name')
parser.add_argument('--rep', choices=['0', '1', '2', '3'], required=True, help = 'Replica')
parser.add_argument('--selection', default = 'Tyr208_His57', help='Distance key defined in distance_defs.py (e.g. Y208-H57)')
parser.add_argument('--list_pairs', default = None, help='List of pairs given by a .pickle file')
parser.add_argument('--out_dir', default='analysis', help = 'Output directory')
parser.add_argument('--top', default='step5_input.psf', help = 'Topology structure')

# parallel controls
parser.add_argument('--chunk_size', type=int, default=20000)
parser.add_argument('--workers', type=int, default=1)

args = parser.parse_args()

root = os.path.abspath(os.path.join(os.getcwd(), '..'))
data_src = os.path.abspath(f'/media/volume/Orai1_rep0/Orai1_analysis')

top = f'system/{args.sys}/{args.top}'
psf = os.path.abspath(os.path.join(data_src, top))

traj_sim = sorted(glob.glob(
    os.path.join(data_src, 'postprocessed_simulation', args.sys, f'Rep-{args.rep}', 'whole/*.xtc')
))

if not traj_sim:
    raise RuntimeError("No trajectories found")

trajs = [os.path.abspath(t) for t in traj_sim]

out_dir = os.path.join(root, args.out_dir, args.sys, 'structural', f'Rep-{args.rep}')
os.makedirs(out_dir, exist_ok=True)

submission_dir = os.path.join(root, 'submission')
jobs_dir = os.path.join(root, 'jobs')

os.makedirs(submission_dir, exist_ok=True)
os.makedirs(jobs_dir, exist_ok=True)

job_name = f'distance_{args.selection}_Rep-{args.rep}'
log_fn = os.path.join(jobs_dir, f'{job_name}.log')
bash_fn = os.path.join(submission_dir, f'{job_name}.sh')

traj_string = " ".join(f'"{t}"' for t in trajs)

out_file = os.path.join(out_dir, f'{args.selection}.h5')

bash_script = f"""#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate openff

echo "===== Distance job started at $(date) ====="
echo "Host: $(hostname)"
echo "CPU cores: $(nproc)"

cd {root}

python -u {root}/scripts/distance.py \\
    --psf "{psf}" \\
    --traj_fns {traj_string} \\
    --selection {args.selection} \\
    --list_pairs {args.list_pairs} \\
    --chunk_size {args.chunk_size} \\
    --workers {args.workers} \\
    --out_fn "{out_file}"

echo "===== Job finished at $(date) ====="
"""


with open(bash_fn, 'w') as f:
    f.write(bash_script)

os.chmod(bash_fn, 0o755)

subprocess.run(
    f'tmux new-session -d -s {job_name} "bash {bash_fn} 2>&1 | tee {log_fn}"',
    shell=True,
    check=True
)


print(f"Job launched in tmux session: {job_name}")
print(f"Log file: {log_fn}")
print(f"Output file: {out_file}")
print(f"Workers: {args.workers}")
print(f"Chunk size: {args.chunk_size}")
print(f"Trajectories found: {len(trajs)}")

