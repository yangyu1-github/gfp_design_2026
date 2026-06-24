# gfp_design_2026
code space for protein design contest 2026
# GFP Design Competition Submission Report

**Team:** ChemEvo-BIT  
**Date:** 2026-06-23  
**Number of sequences:** 6  
**Source campaigns:** Run05 (loop redesign), Run09 (loop-seeded barrel PSSM refinement)  

---

## Design Narrative

Our submission pipeline followed a staged, structure-first strategy:

1. **Loop optimization with RFdiffusion3 and LigandMPNN (Run05).** We treated the GFP beta-barrel as a fixed scaffold and redesigned three surface loops (sfGFP positions 129–147, 188–198, and 209–216) using RFdiffusion3 partial diffusion, followed by LigandMPNN sequence design. Hard constraints preserved the chromophore precursor TYG65–67, the central chromophore helix, catalytic residues, beta-strands, and barrel closure anchors. Candidates were validated with Boltz-2 single-sequence prediction and a TM-align fold audit against PDB 2B3P.

2. **PSSM preparation from loop-screening data.** We used the top Run05 loop designs, together with DMS/alignment data, to compile a position-specific scoring matrix (PSSM) for LigandMPNN. This PSSM guided conservative, evolutionarily informed substitutions at beta-barrel positions surrounding the redesigned loops.

3. **Loop-seeded barrel refinement (Run09).** The two highest-confidence Run05 loop contexts were kept fixed, the low-confidence C-terminal tail after sfGFP coordinate 230 was trimmed, and the remaining beta-barrel positions were redesigned with PSSM-guided LigandMPNN at low temperatures (0.05–0.30). This produced a second-generation set of loop+barrel variants with stronger structural support than direct loop transplantation.

4. **Rosetta FastRelax energy screening for final selection.** All cross-run fold-gated candidates (13 sequences) and selected gate-adjacent controls were subjected to a coordinate-restrained PyRosetta FastRelax ensemble (three decoys per design, 250 minimizer iterations, 0.5 Å C-alpha restraints). The lowest unconstrained `ref2015` energy decoy that stayed within 1.0 Å C-alpha RMSD of the input was selected. We then chose the six designs with the most favorable relaxed energy per residue as our final submission. This Rosetta energy is a structure-conditioned stability proxy, not a literal Cartesian multi-mutation ddG, because the constructs span different loop lengths and mutation backgrounds.

---

## Final Sequences

| Rank | Seq_ID | Source run | Source candidate | Length | Boltz pLDDT | TM2 vs 2B3P | Rosetta REU/residue | Δ REU/residue vs 2B3P |
|------|--------|------------|------------------|--------|-------------|-------------|----------------------|------------------------|
| 1 | chembit_01 | run09 | `run09_candidate_0006_L1_T0p05_score-0.1375_pssm0.0000_lm0.0000_mut29` | 229 | 0.9327 | 0.9195 | -2.9577 | 0.2675 |
| 2 | chembit_02 | run05 | `run05_candidate_0376_run05_three_loop_inputs_run05_three_loop_t3_5_model_4_lm0.4621_mut22` | 239 | 0.8913 | 0.9505 | -2.9292 | 0.2960 |
| 3 | chembit_03 | run05 | `run05_candidate_0152_run05_three_loop_inputs_run05_three_loop_t5_2_model_2_lm0.4833_mut21` | 239 | 0.8962 | 0.9665 | -2.9183 | 0.3069 |
| 4 | chembit_04 | run05 | `run05_candidate_0037_run05_three_loop_inputs_run05_three_loop_t2_3_model_4_lm0.4972_mut20` | 239 | 0.8885 | 0.9649 | -2.9175 | 0.3078 |
| 5 | chembit_05 | run05 | `run05_candidate_0455_run05_three_loop_inputs_run05_three_loop_t3_8_model_1_lm0.4561_mut19` | 239 | 0.8893 | 0.9514 | -2.8767 | 0.3485 |
| 6 | chembit_06 | run05 | `run05_candidate_0196_run05_three_loop_inputs_run05_three_loop_t5_1_model_1_lm0.4789_mut23` | 239 | 0.9035 | 0.9555 | -2.8751 | 0.3501 |

---

## Validation Summary

- **Length:** 229 or 239 amino acids (within 220–250 requirement).
- **Initial residue:** All sequences start with `M`.
- **Alphabet:** Standard 20 amino acids only.
- **Chromophore precursor:** All sequences contain `TYG` at positions 65–67.
- **Exclusion list:** No exact matches to `2026Protein Design/Exclusion_List.csv`.
- **Structural gate:** All selected designs passed the run-specific fold gate (Boltz pLDDT ≥ 0.85 and TM2 ≥ 0.90/0.95).

---

## Selection Rationale

Final selection was driven by the **lowest restrained Rosetta relaxed energy per residue** among candidates that already passed the Boltz structural gate. The top six include:

- One Run09 loop-seeded, coordinate-trimmed control (`run09_candidate_0006`, 229 aa), which combines the high-confidence Run05 L1 loop context with no new barrel mutations and the lowest Rosetta energy in the panel.
- Five Run05 loop-only designs (`run05_candidates 0376, 0152, 0037, 0455, 0196`, 239 aa), which retain the full sfGFP C-terminal tail and show the next-lowest relaxed energies while maintaining high TM2 scores (0.951–0.966).

These six represent the most favorable structure-conditioned energy states discovered in the current campaign.

---

## Caveats

- The Rosetta relaxed energy is an **in-silico stability proxy**, not an experimental thermal stability or brightness measurement.
- Run09 candidates are 229 aa because the C-terminal tail after sfGFP coordinate 230 was intentionally trimmed to avoid low-confidence structure; this is within the competition length window.
- Cartesian multi-mutation ddG could not be completed in the local Rosetta runtime; therefore the energy values are used for relative ranking rather than absolute thermodynamic claims.

---
