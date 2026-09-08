# Stage 1 — Boltz-2 structure prediction

Predicts the Orai1 hexamer bound to the STIM1 CRAC-activation domain (CAD)
in two passes with [LMI4Boltz](https://github.com/tlitfin/lmi4boltz).

## Contents

```
configs/
  orai1_singlechain.yaml                       Boltz-2 input: single Orai1 chain, sequence only
  orai1_hexamer_cad_complex_restrained.yaml     Boltz-2 input: 6x Orai1 + 12x CAD, template + contact restraints
structures/
  orai1_singlechain_model.cif                   Output of the single-chain run
  orai1_singlechain_confidence.json
  orai1_singlechain_template_truncated.cif       Single-chain model, truncated (TM4 extension /
                                                 TM1 cytosolic segment / TM3-TM4 loop removed),
                                                 used as the structural template for the complex run
  orai1_hexamer_cad_complex_model.cif            Output of the restrained complex run
  orai1_hexamer_cad_complex_confidence.json
raw_boltz2_run_output/    (present locally, not tracked in git — see repo-root README)
  boltz_results_orai1.singlechain/
  boltz_results_boltz2-A.trunc/
```

## Method

1. Run Boltz-2 on `configs/orai1_singlechain.yaml` → single-chain Orai1
   model.
2. Manually truncate that model (remove the TM4 extension, the TM1
   cytosolic segment, and the TM3–TM4 loop) → `orai1_singlechain_template_truncated.cif`.
3. Run Boltz-2 on `configs/orai1_hexamer_cad_complex_restrained.yaml`,
   which:
   - defines 6 copies of the Orai1 sequence (chains A–F) and 12 copies of
     the STIM1 CAD sequence (chains G–R, i.e. 2 CAD per Orai1 subunit),
   - supplies the truncated single-chain structure above as a forced
     template for chain A,
   - applies pocket and pairwise-distance restraints encoding known
     Orai1–STIM1 contacts (e.g. Orai1 residue 258 (Tyr) – STIM1 CAD residue
     398 (His); numbering in the config is local to each chain's sequence
     and is annotated in comments with the canonical residue numbers).

Confidence scores for both runs are in the corresponding
`*_confidence.json` file (`confidence_score`, `ptm`/`iptm`, `plddt`, etc.).

(Exact CLI flags depend on the installed [LMI4Boltz](https://github.com/tlitfin/lmi4boltz) — see its own
documentation.)
