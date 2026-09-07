"""
Step 2: Rebuild helix 104-113 as a flexible loop in all Orai1 chains
=====================================================================
Each chain is processed independently through MODELLER LoopModel.
Only the rebuilt loop coordinates are grafted back onto the original
hexamer -- non-loop atoms keep their original positions, preserving
the quaternary structure exactly.

Input:  data/orai1_orai_only.pdb
Output: data/orai1_loop_best.pdb
Next:   python 02b_sequential_loop_minimize.py
"""

import os, shutil, glob

INPUT_PDB      = "data/orai1_model_0.pdb"
BEST_MODEL_OUT = "data/orai1_loop_best.pdb"
LOOP_START     = 105
LOOP_END       = 113
ORAI_CHAINS    = list("ABCDEF")
N_LOOP_MODELS  = 10
MD_LEVEL       = "slow"   # fast | slow | very_slow

# ---------------------------------------------------------------------------
def _fmt_atom_name(name):
    name = name.strip()
    if len(name) == 1: return f" {name}  "
    if len(name) == 2: return f" {name} "
    if len(name) == 3: return f" {name}"
    return name[:4]

# ---------------------------------------------------------------------------
def write_single_chain_pdb(input_pdb, output_pdb, chain_id):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", input_pdb)
    model  = list(struct)[0]
    first_res = last_res = None
    serial = 1
    os.makedirs(os.path.dirname(output_pdb) or ".", exist_ok=True)
    with open(output_pdb, "w") as fh:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                if res.id[0] != " ": continue
                resnum  = res.id[1]
                inscode = res.id[2].strip() or " "
                resname = res.resname[:3]
                cid     = chain.id[0]
                for atom in res:
                    elem   = (atom.element or "").strip()[:2].rjust(2)
                    aname  = _fmt_atom_name(atom.get_name())
                    x, y, z = atom.get_vector().get_array()
                    occ    = atom.occupancy if atom.occupancy is not None else 1.0
                    bfac   = atom.bfactor  if atom.bfactor  is not None else 0.0
                    altloc = atom.get_altloc() or " "
                    line = (
                        f"ATOM  {serial:5d} {aname:<4s}{altloc:1s}"
                        f"{resname:3s} {cid:1s}{resnum:4d}{inscode:1s}   "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}          {elem:>2s}  "
                    )
                    assert len(line) == 80
                    fh.write(line + "\n")
                    serial += 1
                    if first_res is None: first_res = resnum
                    last_res = resnum
            fh.write(f"TER   {serial:5d}      {resname:3s} {cid:1s}{last_res:4d}\n")
            serial += 1
        fh.write("END\n")
    return first_res, last_res

# ---------------------------------------------------------------------------
def write_single_chain_pir(pir_path, pdb_path, chain_id, first_res, last_res):
    from Bio.PDB import PDBParser
    from Bio.SeqUtils import seq1
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", pdb_path)
    res_list = [r for r in list(struct)[0][chain_id] if r.id[0] == " "]
    seq = "".join(seq1(r.resname) for r in res_list) + "*"
    pdb_stem = os.path.splitext(os.path.basename(pdb_path))[0]
    os.makedirs(os.path.dirname(pir_path) or ".", exist_ok=True)
    with open(pir_path, "w") as fh:
        fh.write(">P1;template\n")
        fh.write(f"structureX:{pdb_stem}:{first_res}:{chain_id}:{last_res}:{chain_id}::::\n")
        fh.write(seq + "\n\n")
        fh.write(">P1;target\nsequence:target::::::::\n")
        fh.write(seq + "\n")

# ---------------------------------------------------------------------------
def run_modeller_single_chain(pir_path, pdb_path, chain_id,
                               loop_start, loop_end,
                               n_models, md_level_str, output_dir):
    from modeller import Environ
    from modeller.automodel import LoopModel, refine
    from modeller import Selection

    md_map = {"very_fast": refine.very_fast, "fast": refine.fast,
               "slow": refine.slow, "very_slow": refine.very_slow}
    md_lvl = md_map.get(md_level_str, refine.slow)

    abs_pir        = os.path.abspath(pir_path)
    abs_pdb_dir    = os.path.abspath(os.path.dirname(pdb_path))
    abs_output_dir = os.path.abspath(output_dir)
    orig_dir       = os.path.abspath(".")
    _cid, _ls, _le = chain_id, loop_start, loop_end

    class SingleChainHelixLoop(LoopModel):
        def select_loop_atoms(self):
            for s, e in [(f"{_ls}:{_cid}", f"{_le}:{_cid}"),
                         (f"{_ls}:A",      f"{_le}:A"),
                         (str(_ls),        str(_le))]:
                try:
                    return Selection(self.residue_range(s, e))
                except KeyError:
                    continue
            raise KeyError(f"Cannot find loop {_ls}-{_le}. "
                           f"Chains: {[c.name for c in self.chains]}")

    env = Environ()
    env.io.atom_files_directory = [abs_pdb_dir]
    os.makedirs(abs_output_dir, exist_ok=True)
    os.chdir(abs_output_dir)
    try:
        a = SingleChainHelixLoop(env, alnfile=abs_pir,
                                 knowns="template", sequence="target")
        a.starting_model = a.ending_model = 1
        a.loop.starting_model = 1
        a.loop.ending_model   = n_models
        a.loop.md_level       = md_lvl
        a.make()
    finally:
        os.chdir(orig_dir)

# ---------------------------------------------------------------------------
def fix_chain_id_in_pdb(pdb_path, chain_id):
    lines = []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec in ("ATOM", "TER") and len(line) > 21:
                line = line[:21] + chain_id + line[22:]
            lines.append(line)
    with open(pdb_path, "w") as fh:
        fh.writelines(lines)

# ---------------------------------------------------------------------------
def select_best_loop_model(output_dir, best_out):
    models = glob.glob(os.path.join(output_dir, "target.BL*.pdb"))
    if not models:
        print(f"  ERROR: no loop models found in {output_dir}")
        return False

    log_files = glob.glob(os.path.join(output_dir, "*.log"))
    scores = {}
    for lf in log_files:
        with open(lf) as fh:
            in_summary = header_passed = False
            for line in fh:
                if "Summary of successfully produced models" in line:
                    in_summary = True; header_passed = False; continue
                if in_summary:
                    stripped = line.strip()
                    if not stripped: in_summary = False; continue
                    if "molpdf" in stripped or stripped.startswith("---"):
                        header_passed = True; continue
                    if header_passed:
                        parts = stripped.split()
                        if len(parts) >= 2:
                            try: scores[parts[0]] = float(parts[1])
                            except ValueError: pass

    chain_id = os.path.basename(best_out).split("_")[1]

    if scores:
        best = min(scores, key=scores.get)
        best_path = os.path.join(output_dir, best)
        if not os.path.exists(best_path):
            best_path = models[0]; best = os.path.basename(best_path)
        shutil.copy(best_path, best_out)
        print(f"  Best molpdf={scores.get(best,'?'):.1f} ({best}) -> {best_out}")
    else:
        shutil.copy(models[0], best_out)
        print(f"  (No molpdf scores; copied first model) -> {best_out}")

    fix_chain_id_in_pdb(best_out, chain_id)
    print(f"  Chain ID corrected to '{chain_id}'")
    return True

# ---------------------------------------------------------------------------
def extract_loop_coords(modeller_pdb, loop_start, loop_end):
    """
    Extract atom coordinates for the rebuilt loop region from a MODELLER
    output PDB. Returns {resnum: {atom_name: np.array([x,y,z])}}.
    MODELLER always outputs chain 'A' so we don't filter by chain.
    """
    import numpy as np
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", modeller_pdb)
    model  = list(struct)[0]
    loop_coords = {}
    for chain in model:
        for res in chain:
            if res.id[0] != " ": continue
            resnum = res.id[1]
            if loop_start <= resnum <= loop_end:
                loop_coords[resnum] = {
                    atom.get_name(): atom.get_vector().get_array()
                    for atom in res
                }
    return loop_coords

# ---------------------------------------------------------------------------
def rebuild_hexamer(original_pdb, chain_loop_coords, loop_start, loop_end, output_pdb):
    """
    Graft MODELLER-rebuilt loop coordinates onto the original hexamer.
    Non-loop atoms keep their original positions (hexameric packing preserved).
    Loop atoms (104-113) are replaced with MODELLER's rebuilt coordinates.
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", original_pdb)
    model  = list(struct)[0]
    serial = 1
    n_total = 0
    os.makedirs(os.path.dirname(output_pdb) or ".", exist_ok=True)

    with open(output_pdb, "w") as out:
        for chain in model:
            cid       = chain.id
            loop_data = chain_loop_coords.get(cid, {})
            n_rebuilt = 0
            last_res = None; last_resname = "UNK"

            for res in chain:
                if res.id[0] != " ": continue
                resnum   = res.id[1]
                inscode  = res.id[2].strip() or " "
                resname  = res.resname[:3]
                last_res = resnum; last_resname = resname
                in_loop  = (loop_start <= resnum <= loop_end)
                rebuilt  = loop_data.get(resnum, {}) if in_loop else {}

                atom_order = ["N", "CA", "C", "O", "CB"]
                all_atoms  = {a.get_name(): a for a in res}
                ordered    = [all_atoms[n] for n in atom_order if n in all_atoms]
                rest       = [a for a in res if a.get_name() not in atom_order]

                for atom in ordered + rest:
                    aname  = _fmt_atom_name(atom.get_name())
                    elem   = (atom.element or "").strip()[:2].rjust(2)
                    altloc = atom.get_altloc() or " "
                    occ    = atom.occupancy if atom.occupancy is not None else 1.0
                    bfac   = atom.bfactor  if atom.bfactor  is not None else 0.0
                    raw    = atom.get_name()
                    if in_loop and raw in rebuilt:
                        x, y, z = rebuilt[raw]; n_rebuilt += 1
                    else:
                        x, y, z = atom.get_vector().get_array()
                    line = (
                        f"ATOM  {serial:5d} {aname:<4s}{altloc:1s}"
                        f"{resname:3s} {cid:1s}{resnum:4d}{inscode:1s}   "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}          {elem:>2s}  "
                    )
                    out.write(line + "\n")
                    serial += 1

            if last_res is not None:
                out.write(f"TER   {serial:5d}      {last_resname:3s} {cid:1s}{last_res:4d}\n")
                serial += 1
            n_total += n_rebuilt
            status = f"{n_rebuilt} loop atoms grafted" if loop_data else "original coords kept"
            print(f"  Chain {cid}: {status}")
        out.write("END\n")
    print(f"\n  Hexamer -> {output_pdb}  ({n_total} loop atoms replaced)")

# ---------------------------------------------------------------------------
def verify_loop(pdb_path, loop_start, loop_end):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", pdb_path)
    model = list(s)[0]
    try:    chain = model["A"]
    except: chain = list(model)[0]
    res_dict = {r.id[1]: r for r in chain if r.id[0] == " "}
    print(f"\n  Ca i->i+4 distances in chain A, residues {loop_start}-{loop_end}:")
    any_helix = False
    for i in range(loop_start, loop_end - 3):
        j = i + 4
        if i in res_dict and j in res_dict and "CA" in res_dict[i] and "CA" in res_dict[j]:
            d = (res_dict[i]["CA"].get_vector() - res_dict[j]["CA"].get_vector()).norm()
            tag = "LOOP" if d >= 13 else ("partial" if d >= 10 else "HELIX")
            if d < 10: any_helix = True
            print(f"    {i}->{j}  {d:.2f} A  {tag}")
    if any_helix:
        print("\n  WARNING: still helical. Increase N_LOOP_MODELS or MD_LEVEL='very_slow'.")
    else:
        print("\n  All distances in loop range. Proceed to step 2b.")

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    chain_loop_coords = {}

    for cid in ORAI_CHAINS:
        chain_clean = f"data/chain_{cid}_clean.pdb"
        chain_pir   = f"data/chain_{cid}.pir"
        chain_dir   = f"data/loop_models_chain_{cid}"
        chain_best  = f"data/chain_{cid}_loop_best.pdb"

        print(f"\n{'='*60}\nChain {cid}: writing clean PDB\n{'='*60}")
        first_res, last_res = write_single_chain_pdb(INPUT_PDB, chain_clean, cid)
        print(f"  residues {first_res}-{last_res}")

        print(f"\nChain {cid}: writing PIR")
        write_single_chain_pir(chain_pir, chain_clean, cid, first_res, last_res)
        with open(chain_pir) as f:
            for line in [next(f) for _ in range(2)]:
                print("  ", repr(line.rstrip()))

        print(f"\nChain {cid}: running MODELLER ({N_LOOP_MODELS} models, md={MD_LEVEL})")
        run_modeller_single_chain(chain_pir, chain_clean, cid,
                                  LOOP_START, LOOP_END,
                                  N_LOOP_MODELS, MD_LEVEL, chain_dir)

        print(f"\nChain {cid}: selecting best model")
        ok = select_best_loop_model(chain_dir, chain_best)
        if ok:
            print(f"\nChain {cid}: extracting loop coords ({LOOP_START}-{LOOP_END})")
            loop_coords = extract_loop_coords(chain_best, LOOP_START, LOOP_END)
            chain_loop_coords[cid] = loop_coords
            n_atoms = sum(len(v) for v in loop_coords.values())
            print(f"  {len(loop_coords)} residues, {n_atoms} atoms extracted")
        else:
            print(f"  WARNING: no model for chain {cid} -- original coords kept")

    print(f"\n{'='*60}\nGrafting rebuilt loops onto original hexamer\n{'='*60}")
    rebuild_hexamer(INPUT_PDB, chain_loop_coords, LOOP_START, LOOP_END, BEST_MODEL_OUT)

    print(f"\n{'='*60}\nVerifying loop geometry\n{'='*60}")
    if os.path.exists(BEST_MODEL_OUT):
        verify_loop(BEST_MODEL_OUT, LOOP_START, LOOP_END)

    print(f"\nDone. Output: {BEST_MODEL_OUT}")
    print("Next: python 02b_sequential_loop_minimize.py")
