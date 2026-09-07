"""
Step 2 (assembly): one best loop pose, replicated onto all six Orai1 chains
===========================================================================
MODELLER was run per chain in step 2a, giving six independent rebuilds of the
TM2-TM3 loop.  Those six conformations are NOT equivalent -- they differ by
3-11 A RMSD -- so rather than keeping each chain's own rebuild, this script
takes ONE pose and replicates it onto all six chain positions.  The result is a
C6-symmetric hexamer built from a single, vetted loop conformation.

Pose selection
--------------
SELECTED_CHAIN pins the winner (chain E, chosen after visual inspection:
refined loop, no clashes).  Set it to None to select automatically by RMSD
medoid: every candidate pair is body-superposed on the CA atoms outside the
rebuilt region, the RMSD is measured over that region, and the candidate with
the lowest mean distance to all others -- the most representative pose in the
ensemble -- wins.  The full pairwise matrix is printed either way, so a pinned
choice can always be checked against the medoid.

Replication
-----------
  MODE = "graft_loop"   (default)
      For each chain, the winning model is fitted locally on the backbone of
      the residues flanking GRAFT_START..GRAFT_END, then only that window is
      copied in.  The rest of the hexamer -- including every subunit interface
      and the whole C-terminal CAD-binding region -- keeps its original
      coordinates.
  MODE = "whole_chain"
      The winning model is fitted on all non-loop CA atoms and kept entirely.
      Perfectly symmetric, but the model's body drifts ~2 A from the reference,
      which re-packs the subunit interfaces AND destroys the Orai1-CAD
      interface downstream (measured: 637 Orai1-CAD clashes vs 3).  Kept only
      for comparison.

Why GRAFT_START is 105 and not 110
----------------------------------
MODELLER was asked to rebuild 110-113, but its optimisation also moved the
flanking turns: per-residue backbone deviation from the reference is 3.4-7.2 A
for residues 105-112 and drops to ~1.3 A from 113 onward.  Grafting only
110-113 therefore leaves a 3.8 A gap at the 109-110 junction.  Grafting the
whole displaced stretch 105-113 closes both junctions to normal peptide bonds
(1.37 A and 1.31 A).

Input:  data/orai1_model_0.pdb              (reference frame + chains A-F)
        data/chain_{A..F}_loop_best.pdb     (candidate poses)
Output: data/orai1_refined_loop.pdb
Next:   python 03_rebuild_complex.py
"""

import os
import itertools
import numpy as np

# --------------------------------------------------------------------------- config
REFERENCE_PDB  = "data/orai1_model_0.pdb"
CANDIDATE_PDB  = "data/chain_{cid}_loop_best.pdb"
OUTPUT_PDB     = "data/orai1_refined_loop.pdb"

ORAI_CHAINS    = list("ABCDEF")
CANDIDATES     = list("ABCDEF")      # candidate files to consider

SELECTED_CHAIN = "E"                 # pinned winner; None -> RMSD medoid
SELECT_METRIC  = "window"            # "window" (rebuilt region) or "global" (all CA)

LOOP_START     = 110                 # region MODELLER was asked to rebuild (reporting only)
LOOP_END       = 113
GRAFT_START    = 105                 # region actually transplanted
GRAFT_END      = 113
ANCHOR_SPAN    = 4                   # residues each side used for the local fit
FIT_MARGIN     = 2                   # extra residues each side excluded from the global fit

MODE           = "graft_loop"        # "graft_loop" or "whole_chain"
KEEP_OXT       = True

CLASH_CUTOFF   = 2.0


# --------------------------------------------------------------------------- pdb io
def read_atoms(path, chain_id=None):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            cid = line[21]
            if chain_id is not None and cid != chain_id:
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


def objective_function(path):
    with open(path) as fh:
        for line in fh:
            if "OBJECTIVE FUNCTION" in line:
                try:
                    return float(line.split(":")[-1])
                except ValueError:
                    return None
            if line.startswith("ATOM"):
                break
    return None


def atom_line(serial, a, chain_id):
    return (f"ATOM  {serial:5d} {a['name']:<4s}{a['altloc']:1s}"
            f"{a['resname']:>3s} {chain_id:1s}{a['resseq']:4d}{a['icode']:1s}   "
            f"{a['xyz'][0]:8.3f}{a['xyz'][1]:8.3f}{a['xyz'][2]:8.3f}"
            f"{a['occ']:6.2f}{a['bfac']:6.2f}          {a['element']:>2s}  ")


def index(atoms):
    return {(a["resseq"], a["name"].strip()): a for a in atoms}


# --------------------------------------------------------------------------- math
def kabsch(mobile, target):
    mc, tc = mobile.mean(0), target.mean(0)
    U, _, Vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1.0, 1.0, d]) @ Vt
    return R, tc - mc @ R


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def body_keys(a_idx, b_idx):
    """CA atoms shared by both, outside the transplanted window +/- FIT_MARGIN."""
    lo, hi = GRAFT_START - FIT_MARGIN, GRAFT_END + FIT_MARGIN
    return [k for k in a_idx if k in b_idx and k[1] == "CA" and not (lo <= k[0] <= hi)]


def anchor_keys(a_idx, b_idx):
    """Backbone atoms of the residues flanking the transplanted window."""
    res = (list(range(GRAFT_START - ANCHOR_SPAN, GRAFT_START)) +
           list(range(GRAFT_END + 1, GRAFT_END + 1 + ANCHOR_SPAN)))
    return [(r, n) for r in res for n in ("N", "CA", "C", "O")
            if (r, n) in a_idx and (r, n) in b_idx]


def superpose(mobile_idx, target_idx, keys):
    M = np.array([mobile_idx[k]["xyz"] for k in keys])
    T = np.array([target_idx[k]["xyz"] for k in keys])
    R, t = kabsch(M, T)
    return R, t, rmsd(M @ R + t, T)


# --------------------------------------------------------------------------- pose selection
def pose_distance(a_idx, b_idx, metric="window"):
    """Body-superpose b onto a, then RMSD over the rebuilt window (or all CA)."""
    R, t, _ = superpose(b_idx, a_idx, body_keys(a_idx, b_idx))
    if metric == "global":
        sel = [k for k in a_idx if k in b_idx and k[1] == "CA"]
    else:
        sel = [k for k in a_idx if k in b_idx and GRAFT_START <= k[0] <= GRAFT_END]
    B = np.array([b_idx[k]["xyz"] for k in sel]) @ R + t
    A = np.array([a_idx[k]["xyz"] for k in sel])
    return rmsd(A, B)


def rank_poses(cand_idx, metric):
    n = len(CANDIDATES)
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        d = pose_distance(cand_idx[CANDIDATES[i]], cand_idx[CANDIDATES[j]], metric)
        D[i, j] = D[j, i] = d
    return D, D.sum(1) / (n - 1)


# --------------------------------------------------------------------------- geometry checks
def grid_clashes(coords, labels, cutoff):
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
            cnd = near[(near > i) & (labels[near] != labels[i])]
            if cnd.size:
                d = np.linalg.norm(coords[cnd] - coords[i], axis=1)
                sel = d < cutoff
                hits.extend((i, int(j), float(dd)) for j, dd in zip(cnd[sel], d[sel]))
    return hits


# --------------------------------------------------------------------------- main
def main():
    print("=" * 76)
    print("Step 2: one loop pose, replicated onto all six Orai1 chains")
    print(f"        graft {GRAFT_START}-{GRAFT_END} | anchors +/-{ANCHOR_SPAN} | "
          f"metric {SELECT_METRIC} | mode {MODE}")
    print("=" * 76)

    ref = read_atoms(REFERENCE_PDB)
    ref_by_chain = {}
    for a in ref:
        ref_by_chain.setdefault(a["chain"], []).append(a)
    for c in ORAI_CHAINS:
        if c not in ref_by_chain:
            raise ValueError(f"chain {c} absent from {REFERENCE_PDB}")

    # ---- candidates -------------------------------------------------------
    cand, cand_idx, objf = {}, {}, {}
    for c in CANDIDATES:
        p = CANDIDATE_PDB.format(cid=c)
        if not os.path.exists(p):
            raise FileNotFoundError(p)
        atoms = read_atoms(p)
        if not KEEP_OXT:
            atoms = [a for a in atoms if a["name"].strip() != "OXT"]
        cand[c], cand_idx[c], objf[c] = atoms, index(atoms), objective_function(p)

    print("\n  Candidate poses:")
    for c in CANDIDATES:
        o = "" if objf[c] is None else f"   MODELLER objective {objf[c]:9.2f}"
        print(f"    {CANDIDATE_PDB.format(cid=c):38s} {len(cand[c]):5d} atoms{o}")

    # ---- pairwise RMSD ----------------------------------------------------
    D, mean = rank_poses(cand_idx, SELECT_METRIC)
    print(f"\n  Pairwise RMSD over residues {GRAFT_START}-{GRAFT_END} "
          f"after body superposition (A):")
    print("        " + "  ".join(f"{c:>5s}" for c in CANDIDATES) + "     mean")
    for i, c in enumerate(CANDIDATES):
        print(f"     {c}  " + "  ".join(f"{D[i, j]:5.2f}" for j in range(len(CANDIDATES)))
              + f"    {mean[i]:5.2f}")
    medoid = CANDIDATES[int(np.argmin(mean))]
    print(f"\n  RMSD medoid (lowest mean distance to all others): chain {medoid}")

    winner = SELECTED_CHAIN or medoid
    if winner not in cand_idx:
        raise ValueError(f"SELECTED_CHAIN={winner!r} is not among {CANDIDATES}")
    how = "pinned" if SELECTED_CHAIN else "medoid"
    extra = "" if SELECTED_CHAIN is None else f"; medoid would be {medoid}"
    obj = "" if objf[winner] is None else f", objective {objf[winner]:.2f}"
    print(f"  Winning pose: chain {winner} ({how}{extra})   "
          f"mean RMSD to others {mean[CANDIDATES.index(winner)]:.2f} A{obj}")

    win_idx = cand_idx[winner]

    # ---- replicate onto every chain ---------------------------------------
    print(f"\n  Replicating the chain {winner} pose onto chains {''.join(ORAI_CHAINS)}:")
    built = {}
    for cid in ORAI_CHAINS:
        tgt = [dict(a) for a in ref_by_chain[cid]]
        tgt_idx = index(tgt)

        if MODE == "graft_loop":
            keys = anchor_keys(tgt_idx, win_idx)
            if len(keys) < 8:
                raise ValueError(f"chain {cid}: only {len(keys)} anchor atoms")
            R, t, fit = superpose(win_idx, tgt_idx, keys)
            moved = {k: win_idx[k]["xyz"] @ R + t for k in win_idx
                     if GRAFT_START <= k[0] <= GRAFT_END}
            n = 0
            for a in tgt:
                k = (a["resseq"], a["name"].strip())
                if k in moved:
                    a["xyz"] = moved[k]
                    n += 1
            built[cid] = tgt
            print(f"    chain {cid}: {n:3d} atoms grafted (res {GRAFT_START}-{GRAFT_END})   "
                  f"anchor fit RMSD {fit:4.2f} A ({len(keys)} atoms)")
        else:
            keys = body_keys(tgt_idx, win_idx)
            R, t, fit = superpose(win_idx, tgt_idx, keys)
            chain = [dict(a) for a in cand[winner]]
            for a in chain:
                a["xyz"] = a["xyz"] @ R + t
            built[cid] = chain
            print(f"    chain {cid}: {len(chain):5d} atoms placed   "
                  f"body fit CA RMSD {fit:4.2f} A ({len(keys)} atoms)")

    # ---- write ------------------------------------------------------------
    os.makedirs(os.path.dirname(OUTPUT_PDB) or ".", exist_ok=True)
    serial = 1
    with open(OUTPUT_PDB, "w") as out:
        out.write(f"REMARK   Orai1 hexamer, TM2-TM3 loop rebuilt "
                  f"(MODELLER target {LOOP_START}-{LOOP_END})\n")
        out.write(f"REMARK   single pose from {CANDIDATE_PDB.format(cid=winner)} ({how}), "
                  f"replicated onto all six chains\n")
        out.write(f"REMARK   mode {MODE}, transplanted window {GRAFT_START}-{GRAFT_END}, "
                  f"reference {os.path.basename(REFERENCE_PDB)}\n")
        for cid in ORAI_CHAINS:
            for a in built[cid]:
                out.write(atom_line(serial, a, cid) + "\n")
                serial += 1
            last = built[cid][-1]
            out.write(f"TER   {serial:5d}      {last['resname']:>3s} "
                      f"{cid:1s}{last['resseq']:4d}\n")
            serial += 1
        out.write("END\n")

    total = sum(len(v) for v in built.values())
    print(f"\n  -> {OUTPUT_PDB}   {len(ORAI_CHAINS)} chains, {total} atoms")

    # ---- verification -----------------------------------------------------
    print("\n" + "=" * 76)
    print("Verification")
    print("=" * 76)

    print(f"\n  C6 symmetry: window RMSD of each chain to chain "
          f"{ORAI_CHAINS[0]} after body fit")
    a_idx = index(built[ORAI_CHAINS[0]])
    print("    " + "   ".join(f"{cid} {pose_distance(a_idx, index(built[cid])):.3f} A"
                              for cid in ORAI_CHAINS[1:]))

    print(f"\n  CA i->i+4 across the rebuilt region (helix ~5-6 A, unwound >= 10 A):")
    for cid in ORAI_CHAINS:
        ca = {a["resseq"]: a["xyz"] for a in built[cid] if a["name"].strip() == "CA"}
        parts = []
        for i in range(GRAFT_START - 1, GRAFT_END + 1):
            if i in ca and i + 4 in ca:
                parts.append(f"{i}->{i+4} {np.linalg.norm(ca[i] - ca[i + 4]):5.2f}")
        print(f"    chain {cid}: " + "  ".join(parts))

    print("\n  Backbone continuity (C(i)-N(i+1); a peptide bond is ~1.33 A):")
    worst = 0.0
    for cid in ORAI_CHAINS:
        C = {a["resseq"]: a["xyz"] for a in built[cid] if a["name"].strip() == "C"}
        N = {a["resseq"]: a["xyz"] for a in built[cid] if a["name"].strip() == "N"}
        bad = [(i, float(np.linalg.norm(C[i] - N[i + 1])))
               for i in sorted(C) if i + 1 in N and np.linalg.norm(C[i] - N[i + 1]) > 1.5]
        for cid_i, d in bad:
            worst = max(worst, d)
        if bad:
            print(f"    chain {cid}: " + ", ".join(f"{i}-{i+1} {d:.2f} A" for i, d in bad[:6]))
        else:
            jn = [(GRAFT_START - 1, GRAFT_START), (GRAFT_END, GRAFT_END + 1)]
            txt = ", ".join(f"{i}-{j} {np.linalg.norm(C[i] - N[j]):.2f} A"
                            for i, j in jn if i in C and j in N)
            print(f"    chain {cid}: continuous  (junctions {txt})")

    print(f"\n  Inter-chain heavy-atom contacts < {CLASH_CUTOFF} A:")
    coords, labels, tags, resnums = [], [], [], []
    for ci, cid in enumerate(ORAI_CHAINS):
        for a in built[cid]:
            if a["element"].strip() == "H":
                continue
            coords.append(a["xyz"])
            labels.append(ci)
            tags.append(f"{cid}/{a['resname'].strip()}{a['resseq']}/{a['name'].strip()}")
            resnums.append(a["resseq"])
    hits = grid_clashes(coords, labels, CLASH_CUTOFF)
    if not hits:
        print("    none")
    else:
        print(f"    {len(hits)} pair(s); worst 8:")
        for i, j, d in sorted(hits, key=lambda h: h[2])[:8]:
            print(f"      {tags[i]:26s} -- {tags[j]:26s}  {d:.2f} A")
    n_win = sum(1 for i, j, _ in hits
                if GRAFT_START <= resnums[i] <= GRAFT_END
                or GRAFT_START <= resnums[j] <= GRAFT_END)
    print(f"    involving the transplanted window {GRAFT_START}-{GRAFT_END}: {n_win}")

    print("\nDone. Next: python 03_rebuild_complex.py")


if __name__ == "__main__":
    main()
