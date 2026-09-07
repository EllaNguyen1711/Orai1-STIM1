import os
import glob
import argparse
import subprocess

parser = argparse.ArgumentParser(description="Featurization")
parser.add_argument("--top", type = str, default = 'step5_input.psf', help = 'Topology of MDS, DEFAULT = step5_input.psf')
parser.add_argument("--sysname", type = str, default = 'WT_boltz2', help = 'DEFAULT = popc_0.1mCaCl2_seqfixed_180Apadding')
parser.add_argument("--rep", choices = ['0', '1', '2', '3'], help = 'Replica for TRAJs')
args = parser.parse_args()

root = os.path.abspath(os.path.join(os.getcwd(), '..'))

bash_fn = f'{root}/submission/5-featurization_Rep-{args.rep}.sh'
job_name = f'5-featurization_Rep-{args.rep}'
log_fn = f'{root}/jobs/5-featurization_Rep-{args.rep}.log'

bash_script = f"""#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate openff

cd {root}

python -u {root}/scripts/featurization.py --top {args.top} --sysname {args.sysname} --rep {args.rep}

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

