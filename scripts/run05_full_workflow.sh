#!/usr/bin/env bash
#
# Run05 full three-loop GFP design workflow.
#
# Workflow:
#   1. Generate RFdiffusion3 inputs for combined-loop partial diffusion.
#   2. Run RFdiffusion3 for 200 loop-perturbed backbones.
#   3. Convert CIF to PDB.
#   4. Filter backbones by TM-align vs 2B3P.
#   5. Run LigandMPNN on accepted backbones, designing only loops:
#        129-147, 188-198, 209-216
#   6. Merge, TYG-correct, hard-filter, pathology-filter, and downselect
#      to about 1000 Boltz-2 candidates balanced across backbones.
#   7. Generate Boltz YAMLs, run Boltz-2 empty-MSA prediction, and rank
#      candidates by pLDDT over the designed loop regions.
#
# Usage:
#   bash run05/scripts/run05_full_workflow.sh
#
# Useful overrides:
#   RUN05_DO_RFD3=0 bash run05/scripts/run05_full_workflow.sh
#   RUN05_DO_MPNN=0 bash run05/scripts/run05_full_workflow.sh
#   RUN05_DO_BOLTZ=0 bash run05/scripts/run05_full_workflow.sh
#
set -euo pipefail

SCRIPT_PATH="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
RUN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT_DIR="$(cd "$RUN_DIR/.." && pwd)"
OUT_DIR="$RUN_DIR/output"

source "$HOME/miniconda3/etc/profile.d/conda.sh"

# ---------------------------------------------------------------------------
# Paths and tools
# ---------------------------------------------------------------------------
REF_PDB="${REF_PDB:-$ROOT_DIR/raw/assets/structures/2B3P.pdb}"
EXCLUSION_LIST="${EXCLUSION_LIST:-$ROOT_DIR/raw/current_project/competition_package/Exclusion_List.csv}"
RFD3_INPUT_JSON="${RFD3_INPUT_JSON:-$OUT_DIR/rfd3_inputs/run05_three_loop_inputs.json}"
RFD3_OUT="${RFD3_OUT:-$OUT_DIR/rfd3_output}"
RFD3_PDB_DIR="${RFD3_PDB_DIR:-$OUT_DIR/rfd3_pdbs}"
BACKBONE_FILTER_CSV="${BACKBONE_FILTER_CSV:-$OUT_DIR/backbone_filter.csv}"
ACCEPTED_BACKBONES="${ACCEPTED_BACKBONES:-$OUT_DIR/accepted_backbones.txt}"
MPNN_OUT="${MPNN_OUT:-$OUT_DIR/ligandmpnn}"
RAW_FASTA="${RAW_FASTA:-$OUT_DIR/raw_mpnn_sequences.fasta}"
PROCESSED_FASTA="${PROCESSED_FASTA:-$OUT_DIR/processed_tyg_sequences.fasta}"
PREFILTER_CSV="${PREFILTER_CSV:-$OUT_DIR/sequence_prefilter.csv}"
BOLTZ_FASTA="${BOLTZ_FASTA:-$OUT_DIR/boltz_candidates.fasta}"
BOLTZ_YAML_DIR="${BOLTZ_YAML_DIR:-$OUT_DIR/boltz_yamls}"
BOLTZ_OUT="${BOLTZ_OUT:-$OUT_DIR/boltz_results}"
BOLTZ_ANALYSIS_CSV="${BOLTZ_ANALYSIS_CSV:-$OUT_DIR/boltz_analysis.csv}"
TOP_FASTA="${TOP_FASTA:-$OUT_DIR/top_candidates_for_multisample.fasta}"
LOOP_WINDOWS="${LOOP_WINDOWS:-129-147,188-198,209-216}"
LOOP_PLDDT_ANALYSIS_CSV="${LOOP_PLDDT_ANALYSIS_CSV:-$OUT_DIR/loop_plddt_analysis.csv}"
TOP_LOOP_PLDDT_FASTA="${TOP_LOOP_PLDDT_FASTA:-$OUT_DIR/top_loop_plddt_candidates.fasta}"

LIGANDMPNN_RUNPY="${LIGANDMPNN_RUNPY:-$HOME/Documents/packages/LigandMPNN/run.py}"
LIGANDMPNN_CKPT="${LIGANDMPNN_CKPT:-$HOME/Documents/packages/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt}"
TMALIGN_BIN="${TMALIGN_BIN:-$HOME/Documents/packages/TMalign}"

# ---------------------------------------------------------------------------
# Environments and hardware
# ---------------------------------------------------------------------------
RFD3_ENV="${RFD3_ENV:-rfd3_env}"
LIGANDMPNN_ENV="${LIGANDMPNN_ENV:-ligandmpnn_env}"
BOLTZ_ENV="${BOLTZ_ENV:-boltz_env}"
ANALYSIS_ENV="${ANALYSIS_ENV:-base}"

RFD3_CUDA_VISIBLE_DEVICES="${RFD3_CUDA_VISIBLE_DEVICES:-0}"
MPNN_CUDA_VISIBLE_DEVICES="${MPNN_CUDA_VISIBLE_DEVICES:-1}"
BOLTZ_CUDA_VISIBLE_DEVICES="${BOLTZ_CUDA_VISIBLE_DEVICES:-1}"

# ---------------------------------------------------------------------------
# Scale and thresholds
# ---------------------------------------------------------------------------
RFD3_TOTAL_BACKBONES="${RFD3_TOTAL_BACKBONES:-200}"
RFD3_PARTIAL_T_LIST="${RFD3_PARTIAL_T_LIST:-1 2 3 5}"
RFD3_DIFFUSION_BS="${RFD3_DIFFUSION_BS:-5}"
# If unset, compute batches to preserve approximately RFD3_TOTAL_BACKBONES.
RFD3_BATCHES="${RFD3_BATCHES:-auto}"

BACKBONE_TM_MIN="${BACKBONE_TM_MIN:-0.95}"
BACKBONE_ALIGNED_MIN="${BACKBONE_ALIGNED_MIN:-210}"
MAX_ACCEPTED_BACKBONES="${MAX_ACCEPTED_BACKBONES:-0}"  # 0 = keep all passing

MPNN_SEQS_PER_BACKBONE="${MPNN_SEQS_PER_BACKBONE:-1000}"
MPNN_BATCH_SIZE="${MPNN_BATCH_SIZE:-10}"
MPNN_NUM_BATCHES="${MPNN_NUM_BATCHES:-$(((MPNN_SEQS_PER_BACKBONE + MPNN_BATCH_SIZE - 1) / MPNN_BATCH_SIZE))}"
MPNN_SEED="${MPNN_SEED:-42}"

BOLTZ_TARGET_TOTAL="${BOLTZ_TARGET_TOTAL:-1000}"
TOP_N_FOR_MULTISAMPLE="${TOP_N_FOR_MULTISAMPLE:-24}"

# Stage toggles
RUN05_DO_RFD3="${RUN05_DO_RFD3:-1}"
RUN05_DO_CONVERT="${RUN05_DO_CONVERT:-1}"
RUN05_DO_BACKBONE_FILTER="${RUN05_DO_BACKBONE_FILTER:-1}"
RUN05_DO_MPNN="${RUN05_DO_MPNN:-1}"
RUN05_DO_MERGE="${RUN05_DO_MERGE:-1}"
RUN05_DO_PREFILTER="${RUN05_DO_PREFILTER:-1}"
RUN05_DO_YAMLS="${RUN05_DO_YAMLS:-1}"
RUN05_DO_BOLTZ="${RUN05_DO_BOLTZ:-1}"
RUN05_DO_ANALYZE="${RUN05_DO_ANALYZE:-1}"

mkdir -p "$OUT_DIR" "$RFD3_OUT" "$RFD3_PDB_DIR" "$MPNN_OUT" "$BOLTZ_YAML_DIR" "$BOLTZ_OUT"

activate_env() {
    local env_name="$1"
    set +u
    conda activate "$env_name"
    set -u
}

log_step() {
    echo
    echo "================================================================"
    echo "$1"
    echo "================================================================"
}

echo "Run05 three-loop workflow"
echo "Root:      $ROOT_DIR"
echo "Run dir:   $RUN_DIR"
echo "Output:    $OUT_DIR"
echo "Loops:     129-147, 188-198, 209-216"
echo "RFD3 n:    $RFD3_TOTAL_BACKBONES"
echo "MPNN n/bb: $MPNN_SEQS_PER_BACKBONE"
echo "Boltz n:   $BOLTZ_TARGET_TOTAL"

# ---------------------------------------------------------------------------
# Step 1: Generate RFdiffusion3 input JSON
# ---------------------------------------------------------------------------
log_step "Step 1: Generate RFdiffusion3 inputs"
activate_env "$ANALYSIS_ENV"
python - <<PY
import json
from pathlib import Path

ref_pdb = Path("$REF_PDB")
out_json = Path("$RFD3_INPUT_JSON")
partial_ts = [int(x) for x in "$RFD3_PARTIAL_T_LIST".split()]
if not partial_ts:
    raise SystemExit("No partial_t values configured")

# 2B3P has residues 65-67 as CRO rather than standard chain-A residues.
# Fix everything except the three designed loop modules:
#   129-147, 188-198, 209-216.
fixed = "A2-64,A68-128,A148-187,A199-208,A217-232,B66"

spec = {}
for t in partial_ts:
    spec[f"run05_three_loop_t{t}"] = {
        "input": str(ref_pdb),
        "contig": fixed,
        "select_fixed_atoms": fixed,
        "partial_t": t,
    }

out_json.parent.mkdir(parents=True, exist_ok=True)
out_json.write_text(json.dumps(spec, indent=2) + "\n", encoding="utf-8")
print(f"Wrote {out_json}")
print(f"Cases: {', '.join(spec)}")
print(f"Fixed atoms: {fixed}")
PY

# ---------------------------------------------------------------------------
# Step 2: RFdiffusion3
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_RFD3" == "1" ]]; then
    log_step "Step 2: RFdiffusion3 combined-loop backbone generation"
    activate_env "$RFD3_ENV"
    export CUDA_VISIBLE_DEVICES="$RFD3_CUDA_VISIBLE_DEVICES"

    n_cases=$(python - <<PY
import json
print(len(json.load(open("$RFD3_INPUT_JSON"))))
PY
)
    if [[ "$n_cases" -lt 1 ]]; then
        echo "ERROR: no RFdiffusion3 cases found in $RFD3_INPUT_JSON" >&2
        exit 1
    fi
    rfd3_bs="$RFD3_DIFFUSION_BS"
    if [[ "$rfd3_bs" -lt 1 ]]; then
        rfd3_bs=1
    fi
    if [[ "$RFD3_BATCHES" == "auto" ]]; then
        rfd3_batches=$(((RFD3_TOTAL_BACKBONES + (n_cases * rfd3_bs) - 1) / (n_cases * rfd3_bs)))
    else
        rfd3_batches="$RFD3_BATCHES"
    fi
    echo "RFD3 cases: $n_cases | diffusion_batch_size: $rfd3_bs | batches: $rfd3_batches"
    echo "Approximate RFD3 backbones: $((n_cases * rfd3_bs * rfd3_batches))"

    rfd3 design \
        out_dir="$RFD3_OUT" \
        inputs="$RFD3_INPUT_JSON" \
        n_batches="$rfd3_batches" \
        diffusion_batch_size="$rfd3_bs" \
        skip_existing=True \
        dump_trajectories=False \
        prevalidate_inputs=True \
        2>&1 | tee "$OUT_DIR/rfd3_design.log"
else
    echo "Skipping RFdiffusion3 (RUN05_DO_RFD3=$RUN05_DO_RFD3)"
fi

# ---------------------------------------------------------------------------
# Step 3: Convert CIF to PDB
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_CONVERT" == "1" ]]; then
    log_step "Step 3: Convert RFdiffusion3 CIF outputs to PDB"
    activate_env "$LIGANDMPNN_ENV"
    python "$ROOT_DIR/scripts/convert_cif.py" \
        --input_dir "$RFD3_OUT" \
        --output_dir "$RFD3_PDB_DIR"
else
    echo "Skipping CIF conversion (RUN05_DO_CONVERT=$RUN05_DO_CONVERT)"
fi

# ---------------------------------------------------------------------------
# Step 4: Backbone filter by TM-align
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_BACKBONE_FILTER" == "1" ]]; then
    log_step "Step 4: Backbone filtering by TM-align"
    activate_env "$ANALYSIS_ENV"
    python - <<PY
import csv
import subprocess
from pathlib import Path

pdb_dir = Path("$RFD3_PDB_DIR")
ref_pdb = "$REF_PDB"
tmalign = "$TMALIGN_BIN"
out_csv = Path("$BACKBONE_FILTER_CSV")
accepted_txt = Path("$ACCEPTED_BACKBONES")
tm_min = float("$BACKBONE_TM_MIN")
aligned_min = float("$BACKBONE_ALIGNED_MIN")
max_keep = int("$MAX_ACCEPTED_BACKBONES")

def run_tmalign(pdb: Path):
    result = subprocess.run([tmalign, str(pdb), ref_pdb], capture_output=True, text=True)
    tm1 = tm2 = rmsd = aligned = 0.0
    for line in result.stdout.splitlines():
        if "TM-score=" in line and "normalized by length of Chain_1" in line:
            tm1 = float(line.split("TM-score=")[1].split()[0])
        elif "TM-score=" in line and "normalized by length of Chain_2" in line:
            tm2 = float(line.split("TM-score=")[1].split()[0])
        elif "Aligned length=" in line:
            aligned = float(line.split("Aligned length=")[1].split(",")[0].strip())
            rmsd = float(line.split("RMSD=")[1].split(",")[0].strip())
    return tm1, tm2, rmsd, aligned

rows = []
for pdb in sorted(pdb_dir.glob("*_model_*.pdb")):
    tm1, tm2, rmsd, aligned = run_tmalign(pdb)
    passed = tm2 >= tm_min and aligned >= aligned_min
    rows.append({
        "pdb": str(pdb),
        "name": pdb.stem,
        "tm1": tm1,
        "tm2": tm2,
        "rmsd": rmsd,
        "aligned": aligned,
        "pass": passed,
    })

rows.sort(key=lambda r: (r["pass"], r["tm2"], -r["rmsd"]), reverse=True)
passing = [r for r in rows if r["pass"]]
if max_keep > 0:
    passing = passing[:max_keep]

out_csv.parent.mkdir(parents=True, exist_ok=True)
with out_csv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["pdb", "name", "tm1", "tm2", "rmsd", "aligned", "pass"])
    writer.writeheader()
    writer.writerows(rows)

with accepted_txt.open("w", encoding="utf-8") as handle:
    for row in passing:
        handle.write(row["pdb"] + "\n")

print(f"Backbones scanned: {len(rows)}")
print(f"Backbones passing: {len(passing)}")
print(f"Wrote {out_csv}")
print(f"Wrote {accepted_txt}")
if not passing:
    raise SystemExit("No backbones passed filter")
PY
else
    echo "Skipping backbone filter (RUN05_DO_BACKBONE_FILTER=$RUN05_DO_BACKBONE_FILTER)"
fi

# ---------------------------------------------------------------------------
# Step 5: LigandMPNN, 1000 sequences per accepted backbone
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_MPNN" == "1" ]]; then
    log_step "Step 5: LigandMPNN loop-local sequence design"
    activate_env "$LIGANDMPNN_ENV"
    export CUDA_VISIBLE_DEVICES="$MPNN_CUDA_VISIBLE_DEVICES"

    if [[ ! -s "$ACCEPTED_BACKBONES" ]]; then
        echo "ERROR: accepted backbone list is missing or empty: $ACCEPTED_BACKBONES" >&2
        exit 1
    fi

    # Protein fixed residues. Residues 65-67 are CRO in the mature 2B3P
    # structure and are handled by TYG insertion during sequence postprocess.
    # LigandMPNN expects explicit residue tokens, not range syntax.
    FIXED_RESIDUES="$(python - <<'PY'
ranges = [(2, 64), (68, 128), (148, 187), (199, 208), (217, 232)]
print(" ".join(f"A{i}" for a, b in ranges for i in range(a, b + 1)))
PY
)"
    echo "Fixed residues for LigandMPNN: $FIXED_RESIDUES"
    echo "Design windows: A129-A147 A188-A198 A209-A216"
    echo "LigandMPNN batches: $MPNN_NUM_BATCHES x $MPNN_BATCH_SIZE = $((MPNN_NUM_BATCHES * MPNN_BATCH_SIZE)) per backbone"

    while IFS= read -r pdb; do
        [[ -z "$pdb" ]] && continue
        base="$(basename "$pdb" .pdb)"
        bb_out="$MPNN_OUT/$base"
        mkdir -p "$bb_out"
        echo "LigandMPNN on $base"
        python "$LIGANDMPNN_RUNPY" \
            --seed "$MPNN_SEED" \
            --model_type ligand_mpnn \
            --pdb_path "$pdb" \
            --out_folder "$bb_out" \
            --fixed_residues "$FIXED_RESIDUES" \
            --ligand_mpnn_use_side_chain_context 1 \
            --batch_size "$MPNN_BATCH_SIZE" \
            --number_of_batches "$MPNN_NUM_BATCHES" \
            --checkpoint_ligand_mpnn "$LIGANDMPNN_CKPT" \
            2>&1 | tee -a "$OUT_DIR/ligandmpnn.log"
    done < "$ACCEPTED_BACKBONES"
else
    echo "Skipping LigandMPNN (RUN05_DO_MPNN=$RUN05_DO_MPNN)"
fi

# ---------------------------------------------------------------------------
# Step 6: Merge LigandMPNN FASTAs
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_MERGE" == "1" ]]; then
    log_step "Step 6: Merge LigandMPNN FASTA outputs"
    activate_env "$ANALYSIS_ENV"
    python "$ROOT_DIR/scripts/gfp_merge_mpnn.py" \
        --mpnn_dir "$MPNN_OUT" \
        --out_fasta "$RAW_FASTA"
else
    echo "Skipping merge (RUN05_DO_MERGE=$RUN05_DO_MERGE)"
fi

# ---------------------------------------------------------------------------
# Step 7: TYG correction, hard filters, and downselect to Boltz set
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_PREFILTER" == "1" ]]; then
    log_step "Step 7: TYG correction, hard filters, and Boltz downselect"
    activate_env "$ANALYSIS_ENV"
    python - <<PY
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

raw_fasta = Path("$RAW_FASTA")
processed_fasta = Path("$PROCESSED_FASTA")
prefilter_csv = Path("$PREFILTER_CSV")
boltz_fasta = Path("$BOLTZ_FASTA")
exclusion_csv = Path("$EXCLUSION_LIST")
target_total = int("$BOLTZ_TARGET_TOTAL")

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
TAIL = "GMDELYK"
DESIGN_WINDOWS = [(129, 147), (188, 198), (209, 216)]
ALLOWED = set()
for a, b in DESIGN_WINDOWS:
    ALLOWED.update(range(a, b + 1))

def read_fasta(path):
    records = []
    header = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(re.sub(r"\\s+", "", line))
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records

def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\\n")

def load_exclusion(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        col = "Sequence" if "Sequence" in reader.fieldnames else reader.fieldnames[0]
        return {row[col].strip().upper() for row in reader}

def template_sequence():
    # Template sequence from 2B3P chain A, with TYG restored and sfGFP tail.
    map3 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    seen = set()
    residues = []
    for line in Path("$REF_PDB").read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[21] != "A":
            continue
        resnum = int(line[22:26])
        key = (resnum, line[26])
        if key in seen:
            continue
        seen.add(key)
        if 65 <= resnum <= 67:
            continue
        residues.append((resnum, map3.get(line[17:20].strip(), "X")))
    seq = "M" + "".join(aa for _, aa in sorted(residues))
    if seq[64:67] != "TYG":
        seq = seq[:64] + "TYG" + seq[64:]
    if not seq.endswith(TAIL):
        seq += TAIL
    return seq

def longest_homopolymer(seq):
    best = cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best if seq else 0

def low_complexity(seq, window=12):
    import math
    if len(seq) < window:
        return 0.0
    lows = 0
    total = 0
    for i in range(len(seq) - window + 1):
        w = seq[i:i+window]
        counts = {}
        for c in w:
            counts[c] = counts.get(c, 0) + 1
        h = 0.0
        for v in counts.values():
            p = v / len(w)
            h -= p * math.log(p, 2)
        total += 1
        if h < 2.0:
            lows += 1
    return lows / total if total else 0.0

def lm_conf(header):
    match = re.search(r"overall_confidence=([0-9.]+)", header)
    return float(match.group(1)) if match else 0.0

def backbone_id(header):
    token = header.split()[0]
    return token.split("__seq")[0]

def fix_tyg_and_tail(seq):
    seq = seq.strip().upper()
    if not seq.startswith("M"):
        seq = "M" + seq
    if seq[64:67] != "TYG":
        if len(seq) <= 236:
            seq = seq[:64] + "TYG" + seq[64:]
        else:
            seq = seq[:64] + "TYG" + seq[67:]
    if not seq.endswith(TAIL) and len(seq) < 235:
        seq += TAIL
    return seq

def mutation_positions(ref, seq):
    return [i for i, (a, b) in enumerate(zip(ref, seq), start=1) if a != b]

template = template_sequence()
exclusion = load_exclusion(exclusion_csv)
records = read_fasta(raw_fasta)

processed = []
rows = []
seen = set()
for idx, (header, raw_seq) in enumerate(records, start=1):
    seq = fix_tyg_and_tail(raw_seq)
    processed.append((header, seq))

    reasons = []
    if not (220 <= len(seq) <= 250):
        reasons.append("length")
    if not seq.startswith("M"):
        reasons.append("start_m")
    if set(seq) - AA20:
        reasons.append("noncanonical")
    if seq[64:67] != "TYG":
        reasons.append("missing_tyg")
    if seq in exclusion:
        reasons.append("exclusion")
    if seq in seen:
        reasons.append("duplicate")
    mut_pos = mutation_positions(template, seq)
    outside = [p for p in mut_pos if p not in ALLOWED]
    if outside:
        reasons.append("outside_window:" + ",".join(map(str, outside[:20])))
    if longest_homopolymer(seq) > 4:
        reasons.append("homopolymer")
    if low_complexity(seq) > 0.10:
        reasons.append("low_complexity")

    bb = backbone_id(header)
    conf = lm_conf(header)
    passed = not reasons
    if passed:
        seen.add(seq)
    rows.append({
        "header": header,
        "backbone": bb,
        "sequence": seq,
        "length": len(seq),
        "lm_confidence": conf,
        "mutation_count": len(mut_pos),
        "mutation_positions": ";".join(map(str, mut_pos)),
        "pass": passed,
        "reason": "|".join(reasons) if reasons else "pass",
    })

write_fasta(processed_fasta, processed)

prefilter_csv.parent.mkdir(parents=True, exist_ok=True)
with prefilter_csv.open("w", newline="", encoding="utf-8") as handle:
    fieldnames = ["header", "backbone", "length", "lm_confidence", "mutation_count", "mutation_positions", "pass", "reason", "sequence"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

passing = [r for r in rows if r["pass"]]
groups = defaultdict(list)
for row in passing:
    groups[row["backbone"]].append(row)
for values in groups.values():
    values.sort(key=lambda r: (r["lm_confidence"], r["mutation_count"]), reverse=True)

selected = []
if groups:
    per_backbone = max(1, math.ceil(target_total / len(groups)))
    for bb in sorted(groups):
        selected.extend(groups[bb][:per_backbone])
    selected.sort(key=lambda r: (r["lm_confidence"], r["mutation_count"]), reverse=True)
    selected = selected[:target_total]

boltz_records = []
for i, row in enumerate(selected, start=1):
    header = (
        f"run05_candidate_{i:04d}_{row['backbone']}"
        f"_lm{row['lm_confidence']:.4f}_mut{row['mutation_count']}"
    )
    boltz_records.append((header, row["sequence"]))
write_fasta(boltz_fasta, boltz_records)

print(f"Raw records: {len(records)}")
print(f"Processed records: {len(processed)} -> {processed_fasta}")
print(f"Passing prefilter: {len(passing)}")
print(f"Backbone groups passing: {len(groups)}")
print(f"Boltz selected: {len(boltz_records)} -> {boltz_fasta}")
print(f"Audit CSV: {prefilter_csv}")
if not boltz_records:
    raise SystemExit("No sequences selected for Boltz")
PY
else
    echo "Skipping prefilter (RUN05_DO_PREFILTER=$RUN05_DO_PREFILTER)"
fi

# ---------------------------------------------------------------------------
# Step 8: Generate Boltz YAMLs
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_YAMLS" == "1" ]]; then
    log_step "Step 8: Generate Boltz-2 YAMLs"
    activate_env "$ANALYSIS_ENV"
    rm -rf "$BOLTZ_YAML_DIR"
    mkdir -p "$BOLTZ_YAML_DIR"
    python - <<PY
import re
from pathlib import Path

fasta = Path("$BOLTZ_FASTA")
yaml_dir = Path("$BOLTZ_YAML_DIR")
yaml_dir.mkdir(parents=True, exist_ok=True)

def read_fasta(path):
    records = []
    header = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records

count = 0
for header, seq in read_fasta(fasta):
    name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", header)[:180]
    (yaml_dir / f"{name}.yaml").write_text(
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {seq}\n"
        "      msa: empty\n",
        encoding="utf-8",
    )
    count += 1

print(f"Generated {count} YAMLs in {yaml_dir}")
if count == 0:
    raise SystemExit("No Boltz YAMLs generated")
PY
else
    echo "Skipping Boltz YAML generation (RUN05_DO_YAMLS=$RUN05_DO_YAMLS)"
fi

# ---------------------------------------------------------------------------
# Step 9: Boltz-2 prediction
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_BOLTZ" == "1" ]]; then
    log_step "Step 9: Boltz-2 empty-MSA prediction"
    activate_env "$BOLTZ_ENV"
    export CUDA_VISIBLE_DEVICES="$BOLTZ_CUDA_VISIBLE_DEVICES"

    # Keep BOLTZ_EXTRA_ARGS configurable because Boltz CLI options can vary
    # across installs.
    BOLTZ_EXTRA_ARGS="${BOLTZ_EXTRA_ARGS:-}"
    boltz predict "$BOLTZ_YAML_DIR" \
        --out_dir "$BOLTZ_OUT" \
        $BOLTZ_EXTRA_ARGS \
        2>&1 | tee "$OUT_DIR/boltz.log"
else
    echo "Skipping Boltz (RUN05_DO_BOLTZ=$RUN05_DO_BOLTZ)"
fi

# ---------------------------------------------------------------------------
# Step 10: Analyze Boltz predictions by designed-loop pLDDT
# ---------------------------------------------------------------------------
if [[ "$RUN05_DO_ANALYZE" == "1" ]]; then
    log_step "Step 10: Rank Boltz-2 predictions by designed-loop pLDDT"
    activate_env "$ANALYSIS_ENV"
    python "$SCRIPT_DIR/evaluate_loop_plddt.py" \
        --boltz-out "$BOLTZ_OUT" \
        --fasta "$BOLTZ_FASTA" \
        --out-csv "$LOOP_PLDDT_ANALYSIS_CSV" \
        --top-fasta "$TOP_LOOP_PLDDT_FASTA" \
        --top-n "$TOP_N_FOR_MULTISAMPLE" \
        --loop-windows "$LOOP_WINDOWS"
else
    echo "Skipping analysis (RUN05_DO_ANALYZE=$RUN05_DO_ANALYZE)"
fi

log_step "Run05 workflow complete"
echo "Main outputs:"
echo "  RFD3 input:          $RFD3_INPUT_JSON"
echo "  Backbone filter:     $BACKBONE_FILTER_CSV"
echo "  Accepted backbones:  $ACCEPTED_BACKBONES"
echo "  Raw MPNN FASTA:      $RAW_FASTA"
echo "  Processed FASTA:     $PROCESSED_FASTA"
echo "  Prefilter CSV:       $PREFILTER_CSV"
echo "  Boltz candidates:    $BOLTZ_FASTA"
echo "  Boltz analysis:      $BOLTZ_ANALYSIS_CSV"
echo "  Top candidates:      $TOP_FASTA"
echo "  Loop pLDDT analysis: $LOOP_PLDDT_ANALYSIS_CSV"
echo "  Top loop pLDDT FASTA: $TOP_LOOP_PLDDT_FASTA"
