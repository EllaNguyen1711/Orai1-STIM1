# Stage 3 — Molecular dynamics: system building, simulation, analysis

The minimized complex from Stage 2 (and related structural variants — see
below) is embedded in a POPC bilayer with [CHARMM-GUI](https://charmm-gui.org)
(CHARMM36m force field, CaCl₂-containing solvent), equilibrated and run in
production with OpenMM, and analyzed for structural stability, gating-site
distances, and mutation ΔΔG.

## Contents

```
system/            One folder per CHARMM-GUI-built system (WT_boltz2, WT_refined_loop,
                    boltz2.pentameric). step5_input.{psf,pdb,crd} are the CHARMM-GUI
                    structure/topology files (~200 MB each) — NOT tracked in git.
                    sysinfo.dat (box dimensions) is tracked as a lightweight record.
restraints/        Per-condition restraint definitions (dihe.txt, lipid_pos.txt,
                    prot_pos.txt) for 5 simulation conditions (see below).
toppar/, toppar.str   CHARMM36m / CGenFF force field files (third-party, see repo-root LICENSE).
inputs/            CHARMM-GUI OpenMM input decks: step6.1-6.6 (equilibration), step7 (production).
scripts/           openmm_run.py + omm_*.py (CHARMM-GUI OpenMM driver, unmodified);
                    postprocessing_MD.py, _RMSD.py, _RMSF.py, _Distance.py/distance.py,
                    featurization.py (project analysis code).
submission/        Job wrappers: submit_1 (postprocessing) .. submit_5 (featurization),
                    submit_6/7 (FoldX repair + mutation scan, for DDG estimation).
analysis_notebooks/   DDG_analysis.ipynb, Structural_analysis.RMSF-RMSD.refined.ipynb,
                    Snapshots_selections.ipynb, Y258-F394-H398.ipynb.
```

## Systems and conditions

Three CHARMM-GUI systems are archived locally under `system/`:

- **`WT_boltz2`** — wild-type complex, pre-loop-refinement (Stage 1 output).
- **`WT_refined_loop`** — wild-type complex, post-loop-refinement (Stage 2
  output; the main production system).
- **`boltz2.pentameric`** — an alternative pentameric assembly, built for
  comparison.

`restraints/` covers five conditions in total — the three above plus the
**Y208A** gate mutant and two membrane-only solvation/padding variants
(`popc_0.1mCaCl2_seqfixed_180Apadding`,
`syst_popc-memonly_0.1mMCaCl2_seqfixed[_largepadding]`) — some of which were
built and simulated on a separate compute cluster and are not archived
under `system/` in this copy of the project.

## Simulation protocol

CHARMM-GUI's standard six-stage equilibration (`inputs/step6.1`–`step6.6`,
NVT/NPT with decreasing restraints) followed by unrestrained production
(`inputs/step7_production.inp`, 0.002 ps timestep, PME electrostatics,
force-switch van der Waals, Langevin thermostat at 303.15 K), run with
`scripts/openmm_run.py`.

## Analysis

- **Post-processing** (`submission/submit_1-postprocessing_MD.py` →
  `scripts/postprocessing_MD.py`) — aligns and strides raw trajectories.
- **RMSD** (`submit_2-rmsd.py` → `_RMSD.py`) — whole-system and per-residue.
- **RMSF** (`submit_3-rmsf.py` → `_RMSF.py`).
- **Gate distance** (`submit_4-Y208.H57.py` → `_Distance.py`) — tracks the
  Orai1 Tyr208–STIM1 His57 (canonical numbering: Y258–H398) distance across
  the six Orai1–CAD interfaces, a proxy for channel gating.
- **Featurization** (`submit_5-featurization.py` → `featurization.py`) —
  PyEMMA-based feature extraction for downstream analysis (e.g. Markov
  state modeling).
- **ΔΔG scanning** (`submit_6-foldxrepair.sh`, `submit_7-foldxmutation.sh`)
  — FoldX `RepairPDB` + `BuildModel` on selected snapshots, summarized in
  `analysis_notebooks/DDG_analysis.ipynb`.
- **Notebooks** — `Structural_analysis.RMSF-RMSD.refined.ipynb` (stability),
  `Snapshots_selections.ipynb` (representative-frame selection),
  `Y258-F394-H398.ipynb` (gate-region residue tracking), `DDG_analysis.ipynb`
  (mutation energetics).

## Data availability

The CHARMM-GUI system files (`system/*/step5_input.{psf,pdb,crd}`, raw
trajectories, and the two off-cluster restraint conditions' system files)
are not included in this repository. They can be rebuilt from the Stage 1/2
structures via the [CHARMM-GUI Membrane Builder](https://charmm-gui.org/?doc=input/membrane.bilayer2component),
using the same lipid composition and ion concentration recorded in each
`restraints/<condition>/` set, or made available on request (see repo-root
README).
