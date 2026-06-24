# Run05 Loop Window Plan From 2B3P Secondary Structure

Source: `run05/struture-analysis`

## Numbering Note
2B3P chain A contains protein residues 2-232. The competition sequence has an
N-terminal M at position 1, so residue numbers are intended to match the raw
sfGFP sequence positions after adding M1. In the PDB, residues 65-67 are
represented as the mature chromophore rather than normal ATOM residues.

## Hard-Fixed Regions

Fix these for both RFdiffusion3 and LigandMPNN:

- All beta strands:
  - 12-22
  - 25-36
  - 40-48
  - 92-100
  - 103-115
  - 118-128
  - 148-155
  - 159-171
  - 174-187
  - 199-208
  - 217-227
- Chromophore/central helix region:
  - 57-71
  - especially 65-67 `TYG`
- Catalytic and chromophore-neighbor residues:
  - 96
  - 148
  - 203
  - 222
- N/C barrel closure anchors:
  - 2-22
  - 217-232

## Loop Window Classification

| Window | sfGFP sequence | Class | Rationale |
|--------|----------------|-------|-----------|
| 2-4 | SKG | Avoid | N-terminal barrel closure; too short |
| 9-11 | TGV | Avoid | Near first beta strand; too short |
| 23-24 | NG | Avoid | Strand connector; too short |
| 37-39 | ATN | Avoid | Strand connector; too short |
| 49-56 | TTGKLPVP | Caution | Long enough, but directly precedes central helix 57-71 |
| 64-78 | LTYGVQCFSRYPDHM | Avoid | Contains chromophore precursor and pocket geometry |
| 88-91 | MPEG | Caution | Near catalytic residue 96; short |
| 101-102 | KD | Avoid | Tight turn adjacent to catalytic strand |
| 116-117 | GD | Avoid | Tight beta turn; too short |
| 129-147 | DFKEDGHKLEYNFNS | Primary module | Contiguous non-beta segment between beta strands 118-128 and 148-155; includes short helix 135-138 and loop 139-147 |
| 129-134 | DFKEDG | Fallback | Use only if full 129-147 module is too unstable |
| 139-147 | HKLEYNFNS | Fallback | Use only if full 129-147 module is too unstable; keep residue 148 fixed |
| 156-158 | KQK | Avoid | Too short and adjacent to beta/core |
| 172-173 | VE | Avoid | Too short between beta strands |
| 188-198 | PIGDGPVLLPD | Primary | Longest exposed loop; good novelty window if barrel stays intact |
| 209-216 | SKDPNEKR | Primary | Exposed C-side loop before final beta strand; preserve 217-227 |
| 228-232 | AGITH | Avoid | C-terminal closure anchor; do not perturb |

## Combined-Loop Campaign

Primary run05 should redesign the three useful surface-loop windows together:

- 129-147
- 188-198
- 209-216

This gives a meaningfully novel surface while keeping the GFP fold encoded by
fixed beta strands, the chromophore helix, and the terminal barrel closure.

Fallback split pilots:
- 129-147 alone
- 188-198 alone
- 209-216 alone
- 129-134 or 139-147 only if the full 129-147 module is specifically unstable

## RFdiffusion3 Guidance

Use partial diffusion only:

- pilot `partial_t`: 1, 2, 3, 5
- keep sequence length fixed
- diffuse the three-loop module together
- generate 200 RFdiffusion3 backbones for the combined-loop campaign
- accept a backbone only if the non-diffused barrel aligns cleanly to 2B3P

Do not use high-noise diffusion until the low-noise combined-loop campaign
shows that Boltz-2 preserves the GFP barrel.

## LigandMPNN Guidance

Use `ligand_mpnn` with chromophore context. Fix all residues except the chosen
loop windows and, optionally, explicitly approved exposed flanking residues.

For the main campaign, design residues 129-147, 188-198, and 209-216 only.
Keep adjacent beta-strand anchors fixed:

- 118-128 and 148-155 around window 129-147
- 174-187 and 199-208 around window 188-198
- 199-208 and 217-227 around window 209-216

Generate 1000 LigandMPNN sequences per accepted RFdiffusion3 backbone.

Reject any sequence with changes outside the declared design window.

## Scale And Funnel

1. RFdiffusion3: 200 combined-loop backbones.
2. Backbone filter: barrel alignment, fixed-region RMSD, chromophore/catalytic
   geometry, no loop clashes.
3. LigandMPNN: 1000 sequences per retained backbone.
4. Sequence prefilter: hard contest checks, exact `TYG`, allowed-window-only
   mutations, pathology filters, exclusion-list uniqueness.
5. Boltz-2 first pass: about 1000 total sequences, balanced across retained
   backbones so one backbone cannot dominate.
6. Boltz-2 multi-sample: best 12-24 candidates.
7. Final panel: 6 candidates with novelty, fold confidence, and sequence
   diversity.

## First-Pass Acceptance Gates

For each candidate:

- exact `TYG` at 65-67
- no mutation outside declared loop/design window
- no exclusion-list exact hit
- Boltz-2 empty-MSA single sample:
  - TM-2 >= 0.95 preferred
  - pLDDT >= 0.80 preferred
  - no obvious barrel opening
- multi-sample top candidates:
  - average TM-2 >= 0.90
  - no sample below TM-2 0.75
