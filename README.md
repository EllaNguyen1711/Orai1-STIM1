# Orai1–STIM1 CRAC Channel: Structural Modeling & Molecular Dynamics

Structural modeling and molecular dynamics (MD) simulation pipeline for the
**Orai1 hexameric channel bound to the STIM1 CRAC-activation domain (CAD)**,
the calcium-release-activated calcium (CRAC) channel complex. The pipeline
goes from sequence to a simulation-ready membrane system in three stages:
AI-based structure prediction (Boltz-2), template-guided loop refinement
(MODELLER), and CHARMM-GUI/OpenMM molecular dynamics with downstream
trajectory analysis.

![Boltz-2 model of the Orai1 hexamer bound to twelve STIM1 CAD chains, colored by per-residue pLDDT confidence; top view (left) and side view (right)](docs/figures/orai1_stim1_hexamer_cad_pLDDT.png)

*Boltz-2 prediction of the Orai1 hexamer (chains A–F) bound to twelve copies
of the STIM1 CRAC-activation domain (CAD, chains G–R), colored by
per-residue pLDDT confidence (orange = low, blue = high). Top view (left)
and side view (right). A vector version is available at
[`docs/figures/orai1_stim1_hexamer_cad_pLDDT.pdf`](docs/figures/orai1_stim1_hexamer_cad_pLDDT.pdf).*

## Repository layout

```
.
├── 01_boltz2_structure_prediction/   Stage 1 — Boltz-2 structure prediction
├── 02_modeller_loop_refinement/      Stage 2 — MODELLER loop rebuild & refinement
├── 03_molecular_dynamics/            Stage 3 — CHARMM-GUI system build, OpenMM MD, analysis
└── docs/figures/                     Figures for this README
```

Each stage folder has its own `README.md` with details specific to that
stage. The short version of the pipeline is below.

## Pipeline overview

### Stage 1 — Boltz-2 structure prediction (`01_boltz2_structure_prediction/`)

1. **Single-chain model.** A single Orai1 subunit is predicted from sequence
   alone with [Boltz-2](https://github.com/jwohlwend/boltz)
   (`configs/orai1_singlechain.yaml` → `structures/orai1_singlechain_model.cif`).
2. **Truncated template.** The single-chain model is truncated (TM4
   extension removed, TM1 cytosolic segment removed, TM3–TM4 loop removed)
   to give a minimal, well-supported template
   (`structures/orai1_singlechain_template_truncated.cif`).
3. **Restrained hexamer + CAD complex.** The truncated template is used to
   guide a Boltz-2 prediction of the full complex: six copies of Orai1
   (chains A–F) and twelve copies of the STIM1 CRAC-activation domain (CAD,
   chains G–R, two CAD per Orai1 subunit), with pocket and pairwise-contact
   restraints encoding known Orai1–STIM1 interface residues (e.g. Orai1
   Tyr258–STIM1 His398)
   (`configs/orai1_hexamer_cad_complex_restrained.yaml` →
   `structures/orai1_hexamer_cad_complex_model.cif`).

The full raw Boltz-2 run output (MSA search results, PAE/PDE confidence
matrices, Lightning logs) is kept locally under `raw_boltz2_run_output/` but
is **not** tracked in this repository — see [Data availability](#data-availability-large-files).

### Stage 2 — MODELLER loop refinement (`02_modeller_loop_refinement/`)

The Boltz-2 complex model has an implausible all-helical TM2–TM3 segment
(residues 105–113 in the single-chain numbering). This stage rebuilds that
segment as a flexible loop and reassembles the complex:

| Script | Purpose |
|---|---|
| `00_cif2pdb.py` | Convert the Stage 1 complex `.cif` to clean, per-chain PDBs (Orai1 chains A–F, CAD chains G–R). |
| `01_modeller_loop_rebuild.py` | Run MODELLER `LoopModel` independently on each of the 6 Orai1 chains, rebuilding residues 105–113 as a coil and generating 10 candidate loop conformations per chain. |
| `02_assemble_refined_loop.py` | Select one representative loop pose (by RMSD medoid across chains, or a pinned choice) and graft it onto all six chains, preserving the hexamer's C6 symmetry and the Orai1–CAD interfaces. |
| `03_rebuild_complex.py` | Merge the loop-refined Orai1 hexamer with the untouched CAD chains back into a single complex. |
| `04_minimize_energy.py` | PDBFixer + OpenMM (Amber14/GBn2 implicit solvent) restrained minimization of the rebuilt complex, with a before/after energy comparison against the pre-refinement baseline. |

Result (from `data/energy_report.txt`): the loop graft is relaxable — after
restrained minimization the refined and baseline structures converge to
comparable potential energies (−101,589.7 vs. −101,978.3 kcal/mol; a
difference of +388.5 kcal/mol, essentially noise relative to the >30,000,000
kcal/mol single-point difference before minimization), and inter-chain
steric clashes drop from 21 to 0.

Bulk MODELLER working files (10 raw candidate models plus internal trace
files per chain, ~42 MB) are kept locally under `data/loop_models_chain_*/`
but are not tracked in git; only the best pose per chain
(`data/chain_*_loop_best.pdb`) and the final assembled/minimized structures
are versioned.

### Stage 3 — Molecular dynamics (`03_molecular_dynamics/`)

The refined complex (and related structures — see below) is built into a
membrane system with [CHARMM-GUI](https://charmm-gui.org) (POPC bilayer,
CaCl₂-containing solvent, CHARMM36 force field) and simulated with OpenMM.

- `system/` — one subfolder per CHARMM-GUI-built system (`WT_boltz2`,
  `WT_refined_loop`, `boltz2.pentameric`), including the full CHARMM-GUI
  structure/topology files (`step5_input.psf/.pdb/.crd`).
- `restraints/` — per-condition restraint definitions (protein position,
  lipid position, dihedral) for different systems.
- `toppar/`, `toppar.str` — CHARMM36/CGenFF force field parameter and
  topology files from CHARMM-GUI (third-party; see [License](#license)).
- `inputs/` — CHARMM-GUI OpenMM input decks for six equilibration stages
  (`step6.1`–`step6.6`) and production (`step7`).
- `scripts/` — the CHARMM-GUI-generated OpenMM driver (`openmm_run.py` and
  `omm_*.py` helper modules) plus post-processing and analysis code
  (trajectory alignment, RMSD, RMSF, inter-residue distances, PyEMMA
  featurization).
- `submission/` — job wrappers around the analysis scripts
  (`submit_1`…`submit_5`) and FoldX repair/mutation scans
  (`submit_6`, `submit_7`) for ΔΔG estimation.
- `analysis_notebooks/` — Jupyter notebooks for RMSD/RMSF structural
  analysis, ΔΔG analysis, representative-snapshot selection, and
  gate-residue distance tracking (e.g. Orai1 Y258 / STIM1 F394–H398).

## Data availability (large files)

The CHARMM-GUI-built membrane systems (`03_molecular_dynamics/system/*/step5_input.{psf,pdb,crd}`) are included in this repository.

Two categories of raw pipeline intermediates are still excluded to keep
that number from growing further:

- Raw Boltz-2 run output (`01_boltz2_structure_prediction/raw_boltz2_run_output/`)
- Raw MODELLER candidate models/trace files (`02_modeller_loop_refinement/data/loop_models_chain_*/`)

These can be regenerated from the tracked inputs (Boltz-2 configs, the
refined/minimized PDBs) following the instructions in each stage's
`README.md`. They are also available on request — see [Contact](#contact).

The figures (Figure 6, 7) in this manuscript depict the selected MD snapshot at 430 ns, 
the structure for which can be found at `docs/structures/0430.pdb`. 

## Software

- [Boltz-2](https://github.com/jwohlwend/boltz) — structure prediction
- [MODELLER](https://salilab.org/modeller/) — comparative loop modeling
- [Biopython](https://biopython.org/), [PDBFixer](https://github.com/openmm/pdbfixer), [OpenMM](https://openmm.org/) — structure preparation, minimization, and MD
- [CHARMM-GUI](https://charmm-gui.org) — membrane system building (CHARMM36 force field)
- [MDAnalysis](https://www.mdanalysis.org/) — trajectory analysis (RMSD, RMSF, distances)
- [PyEMMA](http://emma-project.org/) — featurization
- [FoldX](https://foldxsuite.crg.eu/) — ΔΔG mutation scanning
- Python 3.10

## License

Code in this repository is released under the [MIT License](LICENSE). The
force field files under `03_molecular_dynamics/toppar/` are third-party
CHARMM36/CGenFF parameter files distributed by CHARMM-GUI/the MacKerell
lab under their own academic-use terms and are not covered by this
repository's license.

## Contact

Ella — honghanguyen1711@gmail.com
