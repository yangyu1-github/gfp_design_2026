# Run ID: run05

## Status
Completed

## Date
2026-05-12

## Goal
Generate more novel GFP candidates than conservative point mutants while
avoiding the run04 failure mode. Use RFdiffusion3 only for constrained loop
backbone perturbation across selected surface loops, then use LigandMPNN for
loop-local sequence design.

## Strategic Rationale
The contest entry should not rely on run03-style existing mutation
combinations. Those are structurally safe, but likely too close to existing
GFP designs and prior public design space.

Run04 showed the opposite extreme also fails: when beta-strands are allowed to
mutate freely, sequence metrics can look acceptable while Boltz-2 predicts
non-barrel folds. Therefore run05 should search novelty in loop/backbone
surface regions, not in the global beta-barrel encoding.

## Design Principle
Treat the GFP barrel as a fixed fold scaffold and redesign only surface loop
geometry plus loop-facing residues.

Hard-preserve:
- chromophore precursor `TYG` at residues 65-67;
- central helix / chromophore environment, approximately residues 58-71;
- catalytic residues 96 and 222;
- spectral / chromophore-neighbor residues such as 148 and 203 unless a
  specific hypothesis is being tested;
- beta-strand residues and buried core residues from DSSP/SASA analysis;
- N- and C-terminal barrel closure anchors.

Allowed to change:
- solvent-exposed loop residues;
- loop backbone geometry in the selected surface-loop module;
- conservative loop-proximal side chains if Boltz-2 preserves the barrel.

Primary run05 design module:
- 129-147
- 188-198
- 209-216

These three loop regions should be redesigned together to create a more
meaningful novel surface while preserving all beta strands, core residues,
chromophore chemistry, catalytic residues, and barrel closure anchors.

## Proposed Pipeline
Full workflow script:

```bash
bash run05/scripts/run05_full_workflow.sh
```

Memory-safe defaults:
- RFdiffusion3 now uses `RFD3_DIFFUSION_BS=5` and auto-computes enough
  batches to stay near 200 total backbones.
- LigandMPNN remains at `MPNN_BATCH_SIZE=10` and auto-computes enough batches
  to stay near 1000 sequences per accepted backbone.

If GPU memory still fails, use:

```bash
RFD3_DIFFUSION_BS=1 bash run05/scripts/run05_full_workflow.sh
```

1. **Annotate immutable scaffold**
   - Run DSSP/SASA on `raw/assets/structures/2B3P.pdb`.
   - Mark beta-strands, buried residues, chromophore helix, residues 65-67,
     96, 148, 203, 222, and termini as fixed.
   - Define candidate loop windows only from exposed non-strand segments.
   - Current secondary-structure analysis is in `run05/struture-analysis`;
     derived loop priorities are in `run05/loop-window-plan.md`.

2. **RFdiffusion3 loop perturbation**
   - Use partial diffusion with low noise (`partial_t` pilot: 1, 2, 3, 5).
   - Keep sequence length fixed.
   - Diffuse the three-loop module together: 129-147, 188-198, and 209-216.
   - Reject backbones before MPNN if barrel alignment to 2B3P is poor.

3. **LigandMPNN loop-local sequence design**
   - Use `ligand_mpnn`, not `protein_mpnn`.
   - Include CRO or TYG-compatible chromophore context.
   - Fix all scaffold residues; design only residues 129-147, 188-198, and
     209-216 unless a later explicit SASA/core check expands the design set.
   - Generate 1000 sequences per accepted RFdiffusion3 backbone.

4. **Post-process and hard validation**
   - Enforce `M` start, length 220-250, standard AA alphabet.
   - Enforce exact `TYG` at 65-67.
   - Exclusion-list check early and at final output.
   - Reject any sequence with unintended mutations outside allowed windows.

5. **Structure validation**
   - Boltz-2 empty-MSA only.
   - Multi-sample top candidates, at least 3 samples each.
   - Primary gates:
     - TM-2 vs 2B3P >= 0.95 for single-sample triage;
     - average TM-2 >= 0.90 across multi-sample validation;
     - no catastrophic sample below TM-2 0.75;
     - pLDDT preferably >= 0.80.

6. **Sequence triage**
   - ESM-1v or ESM-2 score is secondary only.
   - Do not rank by ESM before passing structural gates.
   - Add pathology filters: homopolymers, low complexity, charge outliers,
     cysteine/proline surprises in loops.

## Pilot Scale
- RFdiffusion3: generate 200 three-loop backbones.
- Backbone filter: retain only backbones with strong barrel preservation and
  intact chromophore/catalytic geometry.
- LigandMPNN: generate 1000 sequences per retained backbone.
- Sequence prefilter: enforce hard constraints, allowed-window mutation limits,
  pathologies, and exclusion-list uniqueness before expensive structure work.
- Boltz-2 triage set: narrow to about 1000 total sequences across all retained
  backbones.
- Multi-sample validate the best 12-24 after first-pass Boltz-2.
- Select 6 only after novelty and hard-filter audit.

## Stop Criteria
Stop the run05 route if:
- best single-sample TM-2 remains below 0.90 after the combined-loop campaign;
- TYG or chromophore-neighbor geometry is not preserved;
- most accepted sequences contain unintended scaffold mutations;
- Boltz-2 multi-sample variance resembles run_native_mpnn instability.

## Expected Outcome
The best run05 designs should be less novel than run04 but far more foldable:
localized loop novelty on a validated GFP scaffold, with no beta-barrel
sequence collapse.

## Result Summary
Run05 achieved the intended middle ground between conservative point mutation
reuse and run04-style scaffold collapse.

Workflow completion:
- RFdiffusion3 generated 200 three-loop backbones.
- Backbone filtering accepted 118/200 backbones.
- LigandMPNN produced 118,118 loop-local sequence records.
- Sequence prefilter retained 87,914 records; failures were mostly duplicate
  sequences, plus 15 homopolymer failures.
- Boltz-2 triaged 1,000 selected sequences with zero failed predictions.
- Post-analysis generated `run05/output/boltz_analysis.csv` and
  `run05/output/top_candidates_for_multisample.fasta`.
- A pooled inspection bundle was created at `run05/top48_boltz_pdbs/`,
  containing converted PDBs, a rank-matched FASTA, a manifest, top-48 analysis
  CSV, and a PyMOL loader.
- Because run05 is a loop-redesign campaign, TM-score should be interpreted as
  a fold-collapse audit rather than the main ranker. A second pooled inspection
  bundle ranked primarily by Boltz pLDDT was created at
  `run05/top48_plddt_boltz_pdbs/`.
- The Run05 evaluator was updated after completion to rank candidates by
  Boltz pLDDT over the designed loop windows only: 129-147, 188-198, and
  209-216. The standalone script is
  `run05/scripts/evaluate_loop_plddt.py`.
- Loop-specific ranking outputs are
  `run05/output/loop_plddt_analysis.csv` and
  `run05/output/top_loop_plddt_candidates.fasta`.

Boltz-2 structural triage:
- Median TM-2 vs 2B3P: 0.9324.
- Best TM-2 vs 2B3P: 0.9665.
- 55/1,000 candidates reached TM-2 >= 0.95.
- 43/1,000 reached TM-2 >= 0.95 and pLDDT >= 0.80.
- 26/1,000 reached TM-2 >= 0.95 and pLDDT >= 0.85.
- 30/1,000 reached TM-2 >= 0.95, pLDDT >= 0.80, and RMSD <= 1.5 A.
- The top-48 FASTA contains 48 valid 239-aa sequences, all starting with `M`,
  using only standard amino acids, and with no exact exclusion-list hits.

Leading candidates:
- Best by designed-loop pLDDT:
  `run05_candidate_0346...t2_0_model_3_lm0.4639_mut28`, loop mean pLDDT
  0.8970, loop minimum pLDDT 0.8286, whole-protein pLDDT 0.9176,
  confidence 0.8990.
- Next loop-pLDDT-ranked candidates:
  `run05_candidate_0205...t1_9_model_3_lm0.4784_mut24`, loop mean pLDDT
  0.8901, and `run05_candidate_0241...t1_9_model_3_lm0.4745_mut24`,
  loop mean pLDDT 0.8896.
- `run05_candidate_0152...t5_2_model_2_lm0.4833_mut21`:
  TM-2 0.9665, pLDDT 0.8962, confidence 0.8268, RMSD 1.24 A.
- `run05_candidate_0037...t2_3_model_4_lm0.4972_mut20`:
  TM-2 0.9649, pLDDT 0.8885, confidence 0.8869, RMSD 1.45 A.
- `run05_candidate_0173...t5_2_model_2_lm0.4818_mut23`:
  TM-2 0.9635, pLDDT 0.8816, confidence 0.8032, RMSD 1.14 A.
- Balanced high-confidence pick:
  `run05_candidate_0196...t5_1_model_1_lm0.4789_mut23`, TM-2 0.9555,
  pLDDT 0.9035, confidence 0.8888, RMSD 1.34 A.

Interpretation:
Run05 is structurally viable. The combined redesign of loops 129-147,
188-198, and 209-216 produced many candidates that preserve the GFP barrel
while carrying roughly 20-25 mutations concentrated in the designed loop
module. The result is more novel than run03 and far more fold-preserving than
run04.

Recommended next step:
Run multi-sample Boltz-2 validation on the pLDDT-ranked records in
`run05/top48_plddt_boltz_pdbs/top48_ranked_by_plddt.fasta`, using TM-score as
a collapse check rather than a major score. Then select a diverse set of 6 by
balancing pLDDT/confidence, backbone diversity, loop mutation diversity, and
chromophore-pocket geometry inspection.
