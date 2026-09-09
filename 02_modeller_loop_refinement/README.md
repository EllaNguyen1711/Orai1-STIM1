# Stage 2 — MODELLER loop refinement

The Stage 1 Boltz-2 complex model predicts residues 105–113 (the TM2–TM3
connection) as an implausible continuous helix in every Orai1 chain. This
stage rebuilds that stretch as a flexible loop with
[MODELLER](https://salilab.org/modeller/), picks one representative
conformation, grafts it onto the hexamer, reassembles the full complex with
the STIM1 CAD chains, and checks the energetic cost of the graft.

## Pipeline

Run in order from this directory:

```bash
python 00_cif2pdb.py
python 01_modeller_loop_rebuild.py
python 02_assemble_refined_loop.py
python 03_rebuild_complex.py
python 04_minimize_energy.py
```

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `00_cif2pdb.py` | `data/boltz2.cif` (Stage 1 complex model) | `data/orai1_model_0.pdb`, `data/orai1_orai_only.pdb` (chains A–F), `data/orai1_cad_only.pdb` (chains G–R), `data/chain_{A..F}_clean.pdb`, `data/chain_{A..F}.pir` |
| 2 | `01_modeller_loop_rebuild.py` | `data/orai1_model_0.pdb` | 10 MODELLER `LoopModel` candidates per chain in `data/loop_models_chain_{A..F}/` (not tracked in git — see repo-root README), best pose copied to `data/chain_{A..F}_loop_best.pdb` |
| 3 | `02_assemble_refined_loop.py` | `data/chain_{A..F}_loop_best.pdb`, `data/orai1_model_0.pdb` | `data/orai1_refined_loop.pdb` — one chosen loop pose grafted onto all six chains, preserving hexamer symmetry |
| 4 | `03_rebuild_complex.py` | `data/orai1_refined_loop.pdb` (Orai1) + `data/orai1_cad_only.pdb` (CAD) | `data/orai1_complex_refined.pdb` |
| 5 | `04_minimize_energy.py` | `data/orai1_complex_refined.pdb`, baseline `data/orai1_model_0.pdb` | `data/orai1_complex_minimized.pdb` (+hydrogens), `data/orai1_complex_minimized_heavy.pdb` (heavy atoms only — this is the structure carried into Stage 3), `data/energy_report.txt` |

### Loop pose selection (step 3)

MODELLER was run independently per chain, so the six rebuilds are not
identical (3–11 Å RMSD apart). Rather than keep six different loop
conformations, `02_assemble_refined_loop.py` selects **one** pose (chain E,
after visual inspection — no clashes) and grafts only the rebuilt window
(residues 105–113, with a superposition margin) onto all six chains,
keeping every other atom — hexamer packing, the Orai1–CAD interfaces — at
its Stage-1 position. Set `SELECTED_CHAIN = None` in the script to select
automatically by RMSD medoid instead of the pinned choice.

### Energetics (step 5)

PDBFixer completes missing atoms/hydrogens at pH 7.4 (no missing residues
are built), then Amber14 + GBn2 implicit solvent is used for a single-point
energy and a restrained minimization (heavy atoms outside the rebuilt
window harmonically restrained).

## Data availability

`data/loop_models_chain_{A..F}/` are kept locally but excluded from git — see the
repo-root README. Re-running step 2 regenerates them.
