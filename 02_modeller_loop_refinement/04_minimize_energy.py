"""
Step 4: energy of the rebuilt complex, before and after restrained minimization
===============================================================================
Answers "how much strain did the loop graft introduce, and how much of it
relaxes away?" with a real force field rather than a restraint score.

Protocol
--------
1. PDBFixer completes the structure: missing heavy atoms, terminal atoms and
   all hydrogens at pH 7.4.  No missing residues are built -- chain breaks and
   termini are left as they are, so nothing is invented.
2. Amber14 + GBn2 implicit solvent, no periodic box.
3. Single-point energy of the structure as built.
4. Minimization with every heavy atom outside FREE_RANGE harmonically
   restrained (RESTRAINT_K), so the loop relaxes while the hexamer packing and
   the Orai1-CAD interface stay where step 2 and step 3 put them.
5. Single-point energy again, plus how far things actually moved.

BASELINE_PDB is run through exactly the same protocol, so the two numbers are
comparable and the difference is attributable to the graft rather than to the
force field, the added hydrogens or the protocol.

A note on the absolute numbers: this is an implicit-solvent calculation on a
membrane protein, with no lipid and no explicit water.  The absolute energies
are not physically meaningful on their own.  The comparison between the two
structures, and the size of the drop on minimization, are what to read.

Input:  data/orai1_complex_refined.pdb   (step 3)
        data/orai1_model_0.pdb           (baseline, pre-graft)
Output: data/orai1_complex_minimized.pdb        (with hydrogens)
        data/orai1_complex_minimized_heavy.pdb  (heavy atoms only, for the pipeline)
        data/energy_report.txt

Runtime is minutes to tens of minutes on a couple of CPU cores; run it in the
background and watch the log.
"""

import os
import sys
import time
import numpy as np

import openmm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer

# --------------------------------------------------------------------------- config
TARGET_PDB    = "data/orai1_complex_refined.pdb"
BASELINE_PDB  = "data/orai1_model_0.pdb"          # None to skip the comparison
OUT_PDB       = "data/orai1_complex_minimized.pdb"
OUT_HEAVY_PDB = "data/orai1_complex_minimized_heavy.pdb"
REPORT_TXT    = "data/energy_report.txt"

FREE_RANGE    = (105, 113)        # residues free to move (the transplanted window)
FREE_CHAINS   = list("ABCDEF")    # ... in these chains
RESTRAINT_K   = 10.0              # kcal/mol/A^2 on every other heavy atom

FORCEFIELD    = ("amber14-all.xml", "implicit/gbn2.xml")
CUTOFF_NM     = 1.2
MAX_ITER      = 500               # 0 = minimise to convergence
TOLERANCE     = 10.0              # kJ/mol/nm
PH            = 7.4
PLATFORM      = "CPU"

KJ_PER_KCAL   = 4.184
CLASH_CUTOFF  = 2.0


def log(msg, fh=None):
    print(msg, flush=True)
    if fh:
        fh.write(msg + "\n")
        fh.flush()


# --------------------------------------------------------------------------- setup
def prepare(path):
    """PDBFixer: complete heavy atoms and add hydrogens. Returns (topology, positions)."""
    fixer = PDBFixer(filename=path)
    fixer.findMissingResidues()
    fixer.missingResidues = {}                 # do not build anything that is not there
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    n_atoms = sum(len(v) for v in fixer.missingAtoms.values())
    n_term = sum(len(v) for v in fixer.missingTerminals.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(PH)
    return fixer.topology, fixer.positions, n_atoms, n_term


def build_system(topology, positions):
    ff = app.ForceField(*FORCEFIELD)
    system = ff.createSystem(topology,
                             nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=CUTOFF_NM * unit.nanometer,
                             constraints=None,
                             rigidWater=False)
    for f in system.getForces():
        f.setForceGroup(0)

    restraint = openmm.CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter(
        "k", RESTRAINT_K * unit.kilocalorie_per_mole / unit.angstrom ** 2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)

    lo, hi = FREE_RANGE
    n_free = 0
    for atom in topology.atoms():
        if atom.element == app.element.hydrogen:
            continue
        res = atom.residue
        free = (res.chain.id in FREE_CHAINS and lo <= int(res.id) <= hi)
        if free:
            n_free += 1
            continue
        restraint.addParticle(atom.index, positions[atom.index].value_in_unit(unit.nanometer))
    restraint.setForceGroup(1)
    system.addForce(restraint)
    return system, restraint.getNumParticles(), n_free


def energies(context):
    """(physical kcal/mol, restraint kcal/mol)"""
    phys = context.getState(getEnergy=True, groups={0}).getPotentialEnergy()
    rest = context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
    return (phys.value_in_unit(unit.kilojoule_per_mole) / KJ_PER_KCAL,
            rest.value_in_unit(unit.kilojoule_per_mole) / KJ_PER_KCAL)


def heavy_coords(topology, positions):
    idx = [a.index for a in topology.atoms() if a.element != app.element.hydrogen]
    P = np.array([positions[i].value_in_unit(unit.angstrom) for i in idx])
    return idx, P


def clash_count(topology, positions, cutoff=CLASH_CUTOFF):
    """Heavy-atom contacts below cutoff between different chains."""
    coords, labels = [], []
    for a in topology.atoms():
        if a.element == app.element.hydrogen:
            continue
        coords.append(positions[a.index].value_in_unit(unit.angstrom))
        labels.append(a.residue.chain.index)
    coords = np.array(coords)
    labels = np.array(labels)
    keys = np.floor(coords / cutoff).astype(int)
    cells = {}
    for i, k in enumerate(map(tuple, keys)):
        cells.setdefault(k, []).append(i)
    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    n, worst = 0, 99.0
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
                n += int(sel.sum())
                if sel.any():
                    worst = min(worst, float(d[sel].min()))
    return n, worst


# --------------------------------------------------------------------------- run one
def run(path, label, fh, minimize=True, write_to=None):
    log(f"\n{'=' * 76}", fh)
    log(f"{label}   ({path})", fh)
    log("=" * 76, fh)

    t0 = time.time()
    topology, positions, n_missing, n_term = prepare(path)
    n_atoms = topology.getNumAtoms()
    n_h = sum(1 for a in topology.atoms() if a.element == app.element.hydrogen)
    log(f"  prepared: {n_atoms} atoms ({n_h} H added, {n_missing} missing heavy atoms "
        f"and {n_term} terminal atoms completed)   [{time.time() - t0:.0f}s]", fh)

    system, n_restrained, n_free = build_system(topology, positions)
    log(f"  restraints: {n_restrained} heavy atoms at "
        f"{RESTRAINT_K} kcal/mol/A^2, {n_free} free "
        f"(res {FREE_RANGE[0]}-{FREE_RANGE[1]} of chains {''.join(FREE_CHAINS)})", fh)

    integrator = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    platform = openmm.Platform.getPlatformByName(PLATFORM)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)

    e0, r0 = energies(context)
    n0, w0 = clash_count(topology, positions)
    log(f"\n  single point", fh)
    log(f"    potential energy   {e0:14,.1f} kcal/mol", fh)
    log(f"    inter-chain contacts < {CLASH_CUTOFF} A: {n0}"
        f"{'' if n0 == 0 else f', closest {w0:.2f} A'}", fh)

    if not minimize:
        return dict(label=label, e0=e0, e1=None, n0=n0, w0=w0)

    t1 = time.time()
    log(f"\n  minimizing (max {MAX_ITER or 'unlimited'} iterations, "
        f"tolerance {TOLERANCE} kJ/mol/nm) ...", fh)
    openmm.LocalEnergyMinimizer.minimize(
        context, TOLERANCE * unit.kilojoule_per_mole / unit.nanometer, MAX_ITER)
    state = context.getState(getPositions=True)
    new_pos = state.getPositions()
    e1, r1 = energies(context)
    log(f"    done in {time.time() - t1:.0f}s", fh)

    idx, P0 = heavy_coords(topology, positions)
    _, P1 = heavy_coords(topology, new_pos)
    lo, hi = FREE_RANGE
    free_mask = np.array([
        (a.residue.chain.id in FREE_CHAINS and lo <= int(a.residue.id) <= hi)
        for a in topology.atoms() if a.element != app.element.hydrogen])
    d = np.linalg.norm(P1 - P0, axis=1)
    n1, w1 = clash_count(topology, new_pos)

    log(f"\n  after restrained minimization", fh)
    log(f"    potential energy   {e1:14,.1f} kcal/mol", fh)
    log(f"    change             {e1 - e0:14,.1f} kcal/mol", fh)
    log(f"    residual restraint {r1:14,.1f} kcal/mol", fh)
    log(f"    heavy-atom displacement: free window "
        f"rmsd {np.sqrt((d[free_mask] ** 2).mean()):.2f} A / max {d[free_mask].max():.2f} A", fh)
    log(f"                             restrained  "
        f"rmsd {np.sqrt((d[~free_mask] ** 2).mean()):.2f} A / max {d[~free_mask].max():.2f} A", fh)
    log(f"    inter-chain contacts < {CLASH_CUTOFF} A: {n0} -> {n1}"
        f"{'' if n1 == 0 else f', closest {w1:.2f} A'}", fh)

    if write_to:
        with open(write_to, "w") as out:
            app.PDBFile.writeFile(topology, new_pos, out, keepIds=True)
        log(f"    -> {write_to}", fh)
        if OUT_HEAVY_PDB:
            with open(write_to) as src, open(OUT_HEAVY_PDB, "w") as dst:
                for line in src:
                    if line.startswith(("ATOM", "HETATM")):
                        el = line[76:78].strip()
                        nm = line[12:16].strip()
                        if el == "H" or (not el and nm.startswith("H")):
                            continue
                    dst.write(line)
            log(f"    -> {OUT_HEAVY_PDB}  (heavy atoms only)", fh)

    return dict(label=label, e0=e0, e1=e1, n0=n0, n1=n1, w0=w0, w1=w1)


# --------------------------------------------------------------------------- main
def main():
    os.makedirs(os.path.dirname(REPORT_TXT) or ".", exist_ok=True)
    with open(REPORT_TXT, "w") as fh:
        log("Step 4: force-field energy of the rebuilt complex", fh)
        log(f"  force field {' + '.join(FORCEFIELD)}, cutoff {CUTOFF_NM} nm, "
            f"platform {PLATFORM}", fh)

        results = [run(TARGET_PDB, "REFINED  (loop 105-113 grafted)", fh,
                       minimize=True, write_to=OUT_PDB)]
        if BASELINE_PDB and os.path.exists(BASELINE_PDB):
            results.append(run(BASELINE_PDB, "BASELINE (pre-graft, orai1_model_0)", fh,
                               minimize=True, write_to=None))

        log(f"\n{'=' * 76}", fh)
        log("Summary  (kcal/mol)", fh)
        log("=" * 76, fh)
        log(f"  {'structure':38s} {'single point':>14s} {'minimized':>14s} {'drop':>12s}", fh)
        for r in results:
            e1 = f"{r['e1']:14,.1f}" if r["e1"] is not None else " " * 14
            dd = f"{r['e1'] - r['e0']:12,.1f}" if r["e1"] is not None else " " * 12
            log(f"  {r['label']:38s} {r['e0']:14,.1f} {e1} {dd}", fh)
        if len(results) == 2:
            a, b = results
            log(f"\n  refined - baseline, single point : {a['e0'] - b['e0']:+,.1f} kcal/mol", fh)
            if a["e1"] is not None and b["e1"] is not None:
                log(f"  refined - baseline, minimized    : {a['e1'] - b['e1']:+,.1f} kcal/mol", fh)
            log("\n  A positive difference is the strain the graft added; if it largely\n"
                "  disappears after minimization, the graft is relaxable and the loop\n"
                "  conformation is viable.", fh)

        log(f"\nReport written to {REPORT_TXT}", fh)


if __name__ == "__main__":
    main()
