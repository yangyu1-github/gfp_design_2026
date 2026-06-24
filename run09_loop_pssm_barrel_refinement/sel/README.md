# Run09 Selected Fold-Gated Candidates

This folder contains the five Run09 candidates with `fold_gate=True` from `top64_boltz_pdbs/manifest.tsv`.

Files:

- `manifest.tsv`: rank, metrics, copied structure paths, and full sequences for the selected candidates.
- `fold_gated_candidates.fasta`: rank-matched FASTA for the same five candidates.
- `pdbs/`: copied rank-prefixed PDB structures for PyMOL/manual review.
- `cifs/`: copied original Boltz model CIFs.
- `load_fold_gated.pml`: PyMOL loader for the copied PDBs.

Selection rule: copied rows where `fold_gate=True`, preserving the original pLDDT ranking.

Note: rank 004 in the top64 package is not copied because it failed the TM-align fold audit; rank 006 is copied because it is fold-gated.
