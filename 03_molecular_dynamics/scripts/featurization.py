import os
import numpy as np
import pyemma.coordinates as coor
import glob
import argparse
import time

parser = argparse.ArgumentParser(description="Featurization")
parser.add_argument("--top", type=str, default='prot.pdb', help='Topology of MDS')
parser.add_argument("--trajs", nargs='+', required=False, help='List of TRAJs')
parser.add_argument("--selection", type=str, default=None, help='Selection of atoms/residues using MDTraj selection language!')
parser.add_argument("--sysname", type=str, default='popc_0.1mCaCl2_seqfixed_180Apadding')
parser.add_argument("--rep", choices=['0', '1', '2', '3'], help='Replica for TRAJs')
args = parser.parse_args()

start = time.time()
psf = os.path.join('system', args.sysname, args.top)

if args.rep:
    trajs = sorted(glob.glob(f'postprocessed_simulation/{args.sysname}/Rep-{args.rep}/protein/*.xtc'))

else:
    trajs = sorted(args.trajs)

if args.selection:
    selection_string = args.selection
else:
    selection_string = ("protein and (chainid 0 to 5) and "
                        "not (resSeq 150 to 180) and "
                        "not (resSeq 240 to 251) and "
                        "not (resSeq 1 to 35)")

feat = coor.featurizer(psf)
feat.select(selection_string)
feat.add_backbone_torsions()
feat.add_sidechain_torsions(which=['chi1', 'chi2'])

data = coor.source(trajs, features=feat)

pca = coor.pca(data=data, var_cutoff=0.9)

os.makedirs(f'analysis/{args.sysname}/structural/Rep-{args.rep}', exist_ok=True)
outfn = f'analysis/{args.sysname}/structural/Rep-{args.rep}/pca.h5'
pca.write_to_hdf5(outfn)

end = time.time()
print(f'Finished within {(end - start)/60} minutes')
