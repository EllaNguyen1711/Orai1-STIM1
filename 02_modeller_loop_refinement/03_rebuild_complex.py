"""
Step 3: rebuild the full Orai1 - CAD complex
============================================
Combines the loop-refined Orai1 hexamer (step 2) with the untouched CAD
subunits, back into a single coordinate file.  Both parts are already in the
frame of data/orai1_model_0.pdb -- step 2 superposed every MODELLER chain onto
its parent chain there, and the CAD chains were never moved -- so this is a
straight merge with no further fitting.

Input:  data/orai1_refined_loop.pdb   chains A-F  (Orai1, TM2-TM3 loop rebuilt)
        data/orai1_cad_only.pdb       chains G-R  (CAD, unchanged)
Output: data/orai1_complex_refined.pdb
"""

import os
import numpy as np

# --------------------------------------------------------------------------- config
ORAI_PDB      = "data/orai1_refined_loop.pdb"
CAD_PDB       = "data/orai1_cad_only.pdb"
REFERENCE_PDB = "data/orai1_model_0.pdb"          # pre-refinement complex, for comparison
OUTPUT_PDB    = "data/orai1_complex_refined.pdb"

ORAI_CHAINS   = list("ABCDEF")
CAD_CHAINS    = list("GHIJKLMNOPQR")

LOOP_START    = 105                               # transplanted window (matches GRAFT_START in step 2)
LOOP_END      = 113
CONTACT_CUT   = 4.0                               # A, heavy-atom interface contact
CLASH_CUT     = 2.0                               # A, heavy-atom clash


# --------------------------------------------------------------------------- pdb io
def read_atoms(path, chains=None):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            cid = line[21]
            if chains is not None and cid not in chains:
                continue
            out.append({
                "name":    line[12:16],
                "altloc":  line[16],
                "resname": line[17:20],
                "chain":   cid,
                "resseq":  int(line[22:26]),
                "icode":   line[26],
                "xyz":     np.array([float(line[30:38]),
                                     float(line[38:46]),
                                     float(line[46:54])]),
                "occ":     float(line[54:60]) if line[54:60].strip() else 1.0,
                "bfac":    float(line[60:66]) if line[60:66].strip() else 0.0,
                "element": line[76:78] if len(line) > 77 else line[12:14],
            })
    return out


def atom_line(serial, a):
    return (f"ATOM  {serial:5d} {a['name']:<4s}{a['altloc']:1s}"
            f"{a['resname']:>3s} {a['chain']:1s}{a['resseq']:4d}{a['icode']:1s}   "
            f"{a['xyz'][0]:8.3f}{a['xyz'][1]:8.3f}{a['xyz'][2]:8.3f}"
            f"{a['occ']:6.2f}{a['bfac']:6.2f}          {a['element']:>2s}  ")


def by_chain(atoms):
    d = {}
    for a in atoms:
        d.setdefault(a["chain"], []).append(a)
    return d


# --------------------------------------------------------------------------- neighbours
def cell_pairs(coords, labels, cutoff):
    """Pairs (i, j, d) with d < cutoff and labels[i] != labels[j]."""
    coords = np.asarray(coords)
    labels = np.asarray(labels)
    keys = np.floor(coords / cutoff).astype(int)
    cells = {}
    for idx, k in enumerate(map(tuple, keys)):
        cells.setdefault(k, []).append(idx)
    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    hits = []
    for k, members in cells.items():
        near = []
        for o in offs:
            near.extend(cells.get((k[0] + o[0], k[1] + o[1], k[2] + o[2]), ()))
        near = np.array(near, dtype=int)
        if near.size == 0:
            continue
        for i in members:
            cand = near[(near > i) & (labels[near] != labels[i])]
            if cand.size:
                d = np.linalg.norm(coords[cand] - coords[i], axis=1)
                sel = d < cutoff
                hits.extend(zip(cand[sel].tolist(), [i] * int(sel.sum()), d[sel].tolist()))
    return hits


def interface_stats(atoms, group_of, cutoff):
    """Contacts between the two groups defined by group_of(atom) -> 0/1/None."""
    coords, labels, idx = [], [], []
    for n, a in enumerate(atoms):
        g = group_of(a)
        if g is None or a["element"].strip() == "H":
            continue
        coords.append(a["xyz"])
        labels.append(g)
        idx.append(n)
    hits = cell_pairs(coords, labels, cutoff)
    return [(idx[i], idx[j], d) for i, j, d in hits]


# --------------------------------------------------------------------------- main
def main():
    print("=" * 70)
    print("Step 3: rebuilding the full Orai1 - CAD complex")
    print("=" * 70)

    for p in (ORAI_PDB, CAD_PDB):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    orai = read_atoms(ORAI_PDB, chains=set(ORAI_CHAINS))
    cad  = read_atoms(CAD_PDB,  chains=set(CAD_CHAINS))

    oc, cc = by_chain(orai), by_chain(cad)
    missing = [c for c in ORAI_CHAINS if c not in oc] + [c for c in CAD_CHAINS if c not in cc]
    if missing:
        raise ValueError(f"missing chains: {missing}")

    print(f"\n  Orai1  {ORAI_PDB}")
    for c in ORAI_CHAINS:
        rs = sorted({a['resseq'] for a in oc[c]})
        print(f"    chain {c}: {len(oc[c]):5d} atoms, res {rs[0]}-{rs[-1]} ({len(rs)} res)")
    print(f"\n  CAD    {CAD_PDB}")
    for c in CAD_CHAINS:
        rs = sorted({a['resseq'] for a in cc[c]})
        print(f"    chain {c}: {len(cc[c]):5d} atoms, res {rs[0]}-{rs[-1]} ({len(rs)} res)")

    # ---- write --------------------------------------------------------------
    os.makedirs(os.path.dirname(OUTPUT_PDB) or ".", exist_ok=True)
    serial = 1
    with open(OUTPUT_PDB, "w") as out:
        out.write("REMARK   Orai1 (A-F) + CAD (G-R) complex\n")
        out.write(f"REMARK   Orai1 TM2-TM3 loop {LOOP_START}-{LOOP_END} rebuilt "
                  f"(step 2, {os.path.basename(ORAI_PDB)})\n")
        out.write(f"REMARK   CAD taken unchanged from {os.path.basename(CAD_PDB)}\n")
        for c in ORAI_CHAINS + CAD_CHAINS:
            atoms = oc.get(c) or cc.get(c)
            for a in atoms:
                out.write(atom_line(serial, a) + "\n")
                serial += 1
            last = atoms[-1]
            out.write(f"TER   {serial:5d}      {last['resname']:>3s} "
                      f"{c:1s}{last['resseq']:4d}\n")
            serial += 1
        out.write("END\n")

    total = len(orai) + len(cad)
    print(f"\n  -> {OUTPUT_PDB}   "
          f"{len(ORAI_CHAINS) + len(CAD_CHAINS)} chains, {total} atoms")

    # ---- verification -------------------------------------------------------
    print("\n" + "=" * 70)
    print("Verification")
    print("=" * 70)

    merged = read_atoms(OUTPUT_PDB)
    print(f"\n  Re-read {len(merged)} atoms from the output "
          f"({'OK' if len(merged) == total else 'MISMATCH'})")

    # CAD must be bit-identical to its source
    src = {(a["chain"], a["resseq"], a["name"]): a["xyz"] for a in cad}
    dev = max(float(np.abs(a["xyz"] - src[(a["chain"], a["resseq"], a["name"])]).max())
              for a in merged if a["chain"] in CAD_CHAINS)
    print(f"  CAD coordinates unchanged: max deviation {dev:.4f} A")

    # Orai1 <-> CAD interface, before and after the loop rebuild
    grp = lambda a: 0 if a["chain"] in ORAI_CHAINS else (1 if a["chain"] in CAD_CHAINS else None)

    ref = read_atoms(REFERENCE_PDB) if os.path.exists(REFERENCE_PDB) else None
    print(f"\n  Orai1-CAD heavy-atom contacts < {CONTACT_CUT} A:")
    new_ct = interface_stats(merged, grp, CONTACT_CUT)
    print(f"    refined complex : {len(new_ct)}")
    if ref is not None:
        ref_ct = interface_stats(ref, grp, CONTACT_CUT)
        print(f"    reference       : {len(ref_ct)}   "
              f"({len(new_ct) - len(ref_ct):+d})")

    print(f"\n  Orai1-CAD clashes < {CLASH_CUT} A:")
    new_cl = interface_stats(merged, grp, CLASH_CUT)
    print(f"    refined complex : {len(new_cl)}")
    if ref is not None:
        ref_cl = interface_stats(ref, grp, CLASH_CUT)
        print(f"    reference       : {len(ref_cl)}   "
              f"({len(new_cl) - len(ref_cl):+d})")
    for i, j, d in sorted(new_cl, key=lambda h: h[2])[:5]:
        ai, aj = merged[i], merged[j]
        print(f"      {ai['chain']}/{ai['resname'].strip()}{ai['resseq']}/{ai['name'].strip():<4s}"
              f" -- {aj['chain']}/{aj['resname'].strip()}{aj['resseq']}/{aj['name'].strip():<4s}  {d:.2f} A")

    # does the rebuilt loop touch CAD at all?
    loop_grp = lambda a: (0 if (a["chain"] in ORAI_CHAINS and LOOP_START <= a["resseq"] <= LOOP_END)
                          else (1 if a["chain"] in CAD_CHAINS else None))
    loop_ct = interface_stats(merged, loop_grp, CONTACT_CUT)
    print(f"\n  Rebuilt loop ({LOOP_START}-{LOOP_END}) to CAD contacts < {CONTACT_CUT} A: {len(loop_ct)}")

    print(f"\nDone. Output: {OUTPUT_PDB}")


if __name__ == "__main__":
    main()
