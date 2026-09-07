"""
Step 1: Convert Boltz-2 CIF output to PDB format
=================================================
Auto-detects Orai1 chains (~251 residues) and CAD chains (~103 residues),
remaps them to standard IDs (A-F for Orai, G-R for CAD), and writes
clean ATOM-only PDB files (exactly 80 chars/line, no SEQRES/REMARK).

Dependencies:
    pip install biopython

Usage:
    python 01_cif_to_pdb.py
"""

import os
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser

CIF_IN = "data/boltz2.cif"
N_ORAI = 6
N_CAD  = 12
ORAI_LEN_RANGE = (230, 270)
CAD_LEN_RANGE  = (80,  125)
ORAI_IDS = list("ABCDEF")
CAD_IDS  = list("GHIJKLMNOPQR")

def _fmt_atom_name(name):
    name = name.strip()
    if len(name) == 1: return f" {name}  "
    if len(name) == 2: return f" {name} "
    if len(name) == 3: return f" {name}"
    return name[:4]

def write_chains_pdb(struct, id_map, output_pdb, orig_chain_order):
    os.makedirs(os.path.dirname(output_pdb) or ".", exist_ok=True)
    serial = 1
    with open(output_pdb, "w") as fh:
        model = list(struct)[0]
        for orig_id in orig_chain_order:
            if orig_id not in id_map:
                print(f"  WARNING: {orig_id} not in id_map -- skipping")
                continue
            new_id = id_map[orig_id]
            try:
                chain = model[orig_id]
            except KeyError:
                print(f"  WARNING: chain {orig_id} not found -- skipping")
                continue
            last_res = None
            last_resname = "UNK"
            for res in chain:
                if res.id[0] != " ":
                    continue
                resnum   = res.id[1]
                inscode  = res.id[2].strip() or " "
                resname  = res.resname[:3]
                last_res = resnum
                last_resname = resname
                atom_order = ["N", "CA", "C", "O", "CB"]
                all_atoms  = {a.get_name().strip(): a for a in res}
                ordered    = [all_atoms[n] for n in atom_order if n in all_atoms]
                rest       = [a for a in res if a.get_name().strip() not in atom_order]
                for atom in ordered + rest:
                    aname  = _fmt_atom_name(atom.get_name())
                    elem   = (atom.element or "").strip()[:2].rjust(2)
                    altloc = atom.get_altloc() or " "
                    occ    = atom.occupancy if atom.occupancy  is not None else 1.0
                    bfac   = atom.bfactor   if atom.bfactor    is not None else 0.0
                    x, y, z = atom.get_vector().get_array()
                    line = (
                        f"ATOM  {serial:5d} {aname:<4s}{altloc:1s}"
                        f"{resname:3s} {new_id:1s}{resnum:4d}{inscode:1s}   "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}          {elem:>2s}  "
                    )
                    assert len(line) == 80, f"ATOM line {len(line)} chars"
                    fh.write(line + "\n")
                    serial += 1
            if last_res is not None:
                fh.write(f"TER   {serial:5d}      {last_resname:3s} {new_id:1s}{last_res:4d}\n")
                serial += 1
        fh.write("END\n")
    new_ids = [id_map[c] for c in orig_chain_order if c in id_map]
    print(f"  Wrote {len(orig_chain_order)} chains {new_ids} -> {output_pdb}")

if __name__ == "__main__":
    print(f"Parsing {CIF_IN} ...")
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("orai1", CIF_IN)
    model = list(struct)[0]

    print("\nChains found in CIF:")
    chain_nres = {}
    for chain in model:
        nres = sum(1 for r in chain if r.id[0] == " ")
        chain_nres[chain.id] = nres
        print(f"  Chain {chain.id:>3s} : {nres} residues")

    orai_orig = sorted([cid for cid, n in chain_nres.items() if ORAI_LEN_RANGE[0] <= n <= ORAI_LEN_RANGE[1]])
    cad_orig  = sorted([cid for cid, n in chain_nres.items() if CAD_LEN_RANGE[0]  <= n <= CAD_LEN_RANGE[1]])
    other     = sorted(set(chain_nres) - set(orai_orig) - set(cad_orig))

    print(f"\nAuto-classified:")
    print(f"  Orai1 ({len(orai_orig)}): {orai_orig}")
    print(f"  CAD   ({len(cad_orig)}): {cad_orig}")
    if other:
        print(f"  Unclassified ({len(other)}): {other}")
    if len(orai_orig) != N_ORAI:
        print(f"  WARNING: expected {N_ORAI} Orai chains, got {len(orai_orig)}")
    if len(cad_orig) != N_CAD:
        print(f"  WARNING: expected {N_CAD} CAD chains, got {len(cad_orig)}")

    id_map = {}
    for i, orig in enumerate(orai_orig): id_map[orig] = ORAI_IDS[i]
    for i, orig in enumerate(cad_orig):  id_map[orig] = CAD_IDS[i]

    print("\nChain ID remapping:")
    for orig, new in id_map.items():
        print(f"  {orig} -> {new}  ({chain_nres[orig]} residues)")

    print()
    write_chains_pdb(struct, id_map, "data/orai1_model_0.pdb",   orai_orig + cad_orig)
    write_chains_pdb(struct, id_map, "data/orai1_orai_only.pdb", orai_orig)
    write_chains_pdb(struct, id_map, "data/orai1_cad_only.pdb",  cad_orig)
    print("\nDone. Next: python 02_modeller_loop_rebuild.py")
