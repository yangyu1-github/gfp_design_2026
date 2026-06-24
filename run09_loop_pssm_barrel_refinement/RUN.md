# Run09: Loop-Seeded Barrel PSSM Refinement

## Status
Completed. Full LigandMPNN, sequence filtering, Boltz validation, analysis, and
top-64 packaging finished by 2026-05-27.

## Date
2026-05-23 to 2026-05-27

## Goal
Use the top Run05 loop designs as fixed loop-seeded GFP backbones trimmed at
sfGFP coordinate 230, then refine beta-barrel residues with Run07-style
DMS/PSSM-guided LigandMPNN.

Run08 showed that direct loop transplantation into Run07 barrels is too brittle.
Run09 keeps the high-confidence Run05 loop context instead, removes the
low-confidence C-terminal tail after residue 230, and redesigns only
PSSM-supported beta-barrel positions around the fixed loop modules.

## Design Scheme
Input loop-seeded backbones:

- Top 2 records from `run05/output/loop_plddt_analysis.csv`.
- PDBs resolved through `run05/top48_plddt_boltz_pdbs/manifest.tsv`.
- Each donor is globally aligned to sfGFP before trimming or extracting loops.
- Trim boundary is sfGFP coordinate 230, not raw donor sequence index 230.
- Prep found that both top Run05 donors lack sfGFP position 171 and carry a
  compensating `G` insertion in the C-terminal tail. Therefore, coordinate-safe
  trimming to sfGFP positions 1-230 gives 229-aa backbones for the current top
  two donors.
- Removed coordinate tail for both current donors: `THGGMDELYK`.
- The raw-index tail after donor position 230 is `HGGMDELYK`; that is not used
  as the default trim rule because it would retain a residue mapped to sfGFP
  position 231.

Fixed sequence motifs:

- exact `TYG65-67`;
- Run05 loop modules, fixed after sfGFP-to-donor mapping:
  - 129-147
  - 188-198
  - 209-216
- Because sfGFP position 171 is absent in the current top two donors, the
  latter two fixed windows map to donor sequence indices 187-197 and 208-215,
  respectively, rather than raw 188-198 and 209-216.

Editable residue policy:

- Use Run07 conservative PSSM positions only if they fall in the 2B3P beta
  strand/barrel ranges from `run05/struture-analysis`.
- Exclude `65-67`, the fixed Run05 loop windows, and positions after 230.
- Expected editable count: 52 residues.
- The 52 sfGFP editable positions are remapped to each trimmed donor's actual
  sequence/PDB tokens before LigandMPNN. In the current top two donors, positions
  after sfGFP 171 are shifted by one residue in the generated PDB token list.

LigandMPNN scale:

- 2 trimmed backbones.
- Temperatures: `0.05`, `0.10`, `0.15`, `0.20`, `0.30`.
- Per backbone-temperature lane: `10,000` sequences.
- Command scale: `--batch_size 100 --number_of_batches 100`.
- Total requested LigandMPNN designs: 100,000.

Boltz scale:

- Select up to 2,000 candidates after sequence filtering.
- Target 200 candidates per backbone-temperature lane before global fill-in.

## Workflow
Full workflow:

```bash
bash run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_barrel_workflow.sh
```

The workflow is resumable by stage toggles. By default it will request 100,000
LigandMPNN sequences, filter to up to 2,000 Boltz candidates, run a 10-structure
5080 smoke batch, then run 50-structure Boltz chunks with 60-second pauses.

Prepare only trimmed inputs and PSSM artifacts:

```bash
RUN09_DO_MPNN=0 \
RUN09_DO_MERGE=0 \
RUN09_DO_PREFILTER=0 \
RUN09_DO_YAMLS=0 \
RUN09_DO_BOLTZ=0 \
RUN09_DO_ANALYZE=0 \
RUN09_DO_PACKAGE=0 \
bash run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_barrel_workflow.sh
```

Run filtering and downstream stages after LigandMPNN:

```bash
RUN09_DO_PREP=0 \
RUN09_DO_MPNN=0 \
bash run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_barrel_workflow.sh
```

## Filtering And Ranking
Hard sequence filters:

- length exactly matches the prepared trimmed backbone for that donor
  (229 aa for the current top two coordinate-trimmed inputs);
- starts with `M`;
- standard amino acids only;
- exact `TYG65-67`;
- no exact exclusion-list hit;
- no duplicate sequence;
- fixed Run05 loop windows remain exactly unchanged;
- no mutation relative to the prepared Run05 donor outside the declared
  beta-barrel PSSM editable set.

Mutation accounting:

- `mutation_count` is reported versus sfGFP coordinates and includes inherited
  Run05 loop/source differences plus new MPNN edits.
- `new_barrel_mutation_count` counts only edits introduced by LigandMPNN at
  the mapped beta-barrel PSSM sites.
- Loop pLDDT diagnostics use mapped donor sequence indices, not raw sfGFP
  indices, after candidate metadata are available.

Pathology filters:

- homopolymer >4;
- dipeptide repeat >3;
- low-complexity fraction >0.05;
- net charge outside `[-20, 20]`;
- cysteine count >3.

Soft Boltz downselect ranking:

- LigandMPNN confidence;
- ligand confidence;
- PSSM mean;
- closeness to mutation target 40;
- lower new beta-barrel mutation burden;
- sequence diversity.

## Boltz Policy
Boltz-2 validation is empty-MSA only.

- RTX 5080 first:
  - `CUDA_DEVICE_ORDER=PCI_BUS_ID`
  - `CUDA_VISIBLE_DEVICES=1`
  - expected GPU name contains `5080`
- A4000 fallback:
  - `CUDA_VISIBLE_DEVICES=0`
  - expected GPU name contains `A4000`
- Smoke batch: 10 structures.
- Full chunks: 50 structures.
- Pause: 60 seconds between chunks.
- Full and fallback Boltz loops skip chunks whose expected model files already
  exist, so fallback is used for missing or failed chunks rather than blindly
  restarting completed chunks.

Fold gates reported after Boltz:

- pLDDT >= 0.85;
- pLDDT >= 0.90;
- TM2 audit >= 0.75;
- full fold gate: pLDDT >= 0.85, TM2 >= 0.75, aligned residues >= 210,
  RMSD <= 3.5 A.

## Validation Performed
- `python -m py_compile run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_utils.py`
- `bash -n run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_barrel_workflow.sh`
- Prep-only workflow dry run completed.
- Prepared 2 trimmed backbones, both length 229 after sfGFP-coordinate trim.
- Confirmed 52 mapped editable PSSM residues per backbone.
- Ran a two-sequence prefilter smoke test using the prepared trimmed backbone
  sequences; both passed and produced Boltz candidate records.
- Ran a Boltz YAML chunk-generation smoke test with the prepared trimmed
  backbone FASTA; it generated two YAMLs, two one-record chunks, and one smoke
  YAML.

## Execution Summary

LigandMPNN and filtering:

- Raw LigandMPNN records: 100,010.
- Processed records: 100,010.
- Passing sequence-filter records: 12,061.
- Boltz candidates selected: 2,000.
- Lane balancing: 200 selected candidates for each of the 10
  backbone-temperature lanes.
- Editable residues per backbone: 52.
- Main sequence-filter failure modes:
  - MPNN confidence below 0.43: 66,894 records.
  - Duplicate sequence: 52,826 records.
  - Cysteine count >3: 214 records.

Pass rates by loop donor:

- L1 / `run05_candidate_0346...mut28`: 9,374 passing records.
- L2 / `run05_candidate_0205...mut24`: 2,687 passing records.
- Interpretation: the L1 loop context sampled more permissively under the
  Run09 PSSM-barrel policy.

Boltz validation:

- Boltz candidates analyzed: 2,000/2,000.
- pLDDT median: 0.5831.
- pLDDT range: 0.3691-0.9327.
- confidence median: 0.5288.
- TM2 median: 0.6322.
- TM2 range: 0.2922-0.9256.
- loop pLDDT mean median: 0.4762.

Gate counts:

- pLDDT >= 0.85: 12/2,000.
- pLDDT >= 0.90: 1/2,000.
- TM2 >= 0.75: 170/2,000.
- TM2 >= 0.90: 19/2,000.
- pLDDT >= 0.85 and TM2 >= 0.90: 3/2,000.
- pLDDT >= 0.90 and TM2 >= 0.90: 1/2,000.
- Full fold gate: 5/2,000.

Packaged outputs:

- `output/boltz_analysis.csv`
- `output/top_candidates_by_plddt.fasta`
- `output/top_fold_gated_candidates.fasta`
- `top64_boltz_pdbs/manifest.tsv`
- `top64_boltz_pdbs/top_ranked_by_plddt.fasta`
- `top64_boltz_pdbs/load_top_by_plddt.pml`
- `top64_boltz_pdbs/pdbs/`

## Top Candidates

Fold-gated candidates:

| Rank | Candidate | Donor | Temp | pLDDT | Conf | TM2 | RMSD | Aligned | Mutations | New barrel mutations | Loop pLDDT |
|------|-----------|-------|------|-------|------|-----|------|---------|-----------|----------------------|------------|
| 1 | `run09_candidate_0006_L1_T0p05_score-0.1375_pssm0.0000_lm0.0000_mut29` | L1 | 0.05 | 0.9327 | 0.8091 | 0.9195 | 1.88 | 223 | 29 | 0 | 0.9057 |
| 2 | `run09_candidate_0180_L1_T0p05_score2.9568_pssm-0.1410_lm0.4432_mut51` | L1 | 0.05 | 0.8887 | 0.7713 | 0.8851 | 1.97 | 222 | 51 | 23 | 0.8194 |
| 3 | `run09_candidate_1004_L2_T0p05_score-0.1875_pssm0.0000_lm0.0000_mut25` | L2 | 0.05 | 0.8819 | 0.7712 | 0.9256 | 1.97 | 224 | 25 | 0 | 0.8236 |
| 4 | `run09_candidate_0451_L1_T0p15_score2.9874_pssm-0.1544_lm0.4417_mut49` | L1 | 0.15 | 0.8760 | 0.7687 | 0.9114 | 1.67 | 221 | 49 | 21 | 0.8315 |
| 5 | `run09_candidate_0245_L1_T0p10_score2.9822_pssm-0.1744_lm0.4454_mut49` | L1 | 0.10 | 0.8685 | 0.7496 | 0.8698 | 2.07 | 220 | 49 | 21 | 0.8204 |

Notes:

- Candidate `0006` is the strongest structural result, but it has zero new
  barrel mutations and appears to be the trimmed L1 Run05 donor context rather
  than a newly refined barrel design. It is important as a validation/control
  for the coordinate-trimmed Run05 loop context, but it is not the main novelty
  product of Run09.
- Candidate `1004` similarly validates the trimmed L2 donor context with zero
  new barrel mutations.
- The strongest true PSSM-refined candidates are `0180`, `0451`, and `0245`,
  all from the L1 donor and all carrying 21-23 new barrel mutations.
- Low temperature performed best for fold preservation: full fold-gated
  candidates came from temperatures 0.05, 0.10, and 0.15; none came from 0.20
  or 0.30.

## Interpretation

Run09 partially supports the hypothesis.

Positive findings:

- Keeping the Run05 loop context fixed and trimming the low-confidence
  C-terminal tail is structurally viable.
- The coordinate-trimmed top Run05 donors themselves validate well under
  Boltz, especially L1.
- Conservative PSSM-guided barrel refinement can produce fold-gated candidates
  with pLDDT >0.86 and TM2 close to or above 0.90.
- Run09 is a clear improvement over Run08 direct grafting, which produced no
  TM2 >=0.90 grafts.

Limitations and failure modes:

- Most designed candidates still collapsed or lost confident barrel alignment:
  median TM2 was only 0.6322 and median pLDDT was 0.5831.
- The highest pLDDT candidate is a zero-new-barrel-mutation control-like
  sequence, so the best result does not prove that broad beta-barrel editing is
  generally safe.
- Only 5/2,000 candidates passed the full fold gate.
- L2 was much less permissive during sequence filtering and contributed only
  one fold-gated candidate, also with zero new barrel mutations.
- MPNN confidence and duplicate generation were the dominant pre-Boltz filter
  losses.

## Decision

- Keep the top Run09 fold-gated candidates for manual PyMOL review,
  especially `0180`, `0451`, and `0245` as true PSSM-refined L1 designs.
- Treat `0006` and `1004` as coordinate-trimmed Run05 donor controls rather
  than final novelty candidates.
- Do not broaden the editable beta-barrel set from here. The next refinement
  should become more conservative, favoring L1 and a smaller subset of
  high-tolerance barrel positions that appear in successful Run09 candidates.
- Use temperature 0.05-0.15 for the next sequence-design round.
- Before synthesis decisions, inspect the chromophore pocket and barrel
  register manually for the 5 fold-gated candidates.

## Implemented Artifacts
- `run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_barrel_workflow.sh`
- `run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_utils.py`

Prep outputs already generated:

- `input/trimmed_backbones.tsv`
- `input/mapping_qc.csv`
- `input/trimmed_backbones/*.pdb`
- `pssm/barrel_pssm_positions.csv`
- per-backbone LigandMPNN residue/bias/omit JSON files

Expected outputs:

- `input/trimmed_backbones.tsv`
- `input/trimmed_backbones/*.pdb`
- `input/mapping_qc.csv`
- `pssm/redesigned_residues_barrel_pssm.txt`
- `pssm/barrel_pssm_positions.csv`
- `output/raw_mpnn_sequences.fasta`
- `output/sequence_prefilter.csv`
- `output/sequence_prefilter_summary.md`
- `output/boltz_candidates.fasta`
- `output/boltz_candidates_manifest.csv`
- `output/boltz_analysis.csv`
- `output/top_candidates_by_plddt.fasta`
- `output/top_fold_gated_candidates.fasta`
- `top64_boltz_pdbs/`

## Test Plan
Before LigandMPNN:

- Confirm exactly 2 loop donors selected by Run05 loop pLDDT rank.
- Confirm both donor PDBs are found.
- Confirm sequence-to-PDB mapping is unambiguous.
- Confirm trimmed PDB and sequence lengths match the coordinate-aware trim
  result; current top two donors are 229 aa after trimming to sfGFP coordinate
  230.
- Confirm `TYG65-67` and all fixed loop windows map correctly after trimming.
- Confirm editable residue count is 52.

After LigandMPNN:

- Confirm requested raw scale is 100,000.
- Report pass/fail counts by backbone and temperature.
- Report duplicate rate, mutation counts, MPNN confidence, and pathology
  failures.
- Confirm all passing candidates match their prepared donor length and preserve
  fixed loops.

Before Boltz:

- Confirm up to 2,000 selected candidates.
- Confirm lane balancing target of 200 per backbone-temperature lane where
  enough candidates pass.
- Run 5080 preflight and smoke batch.

After Boltz:

- Report pLDDT, confidence, pTM, TM2 audit, RMSD audit, and aligned residues.
- Count pLDDT >= 0.85, pLDDT >= 0.90, TM2 >= 0.75, and full fold-gated
  candidates.
- Package top structures and rank-matched FASTA.
