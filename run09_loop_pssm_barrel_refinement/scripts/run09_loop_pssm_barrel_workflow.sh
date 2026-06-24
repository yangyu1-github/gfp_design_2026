#!/usr/bin/env bash
#
# Run09 loop-seeded barrel PSSM refinement workflow.
#
# Uses the top 2 Run05 loop-pLDDT designs as sfGFP-coordinate 1-230 trimmed
# backbones, then performs Run07-style DMS/PSSM-guided LigandMPNN on
# conservative beta-barrel positions only. Boltz validation is empty-MSA, RTX
# 5080 first, with A4000 fallback.
#
# Usage:
#   bash run09_loop_pssm_barrel_refinement/scripts/run09_loop_pssm_barrel_workflow.sh
#
set -euo pipefail

SCRIPT_PATH="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
RUN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT_DIR="$(cd "$RUN_DIR/.." && pwd)"
INPUT_DIR="${INPUT_DIR:-$RUN_DIR/input}"
PSSM_DIR="${PSSM_DIR:-$RUN_DIR/pssm}"
OUT_DIR="${OUT_DIR:-$RUN_DIR/output}"

source "$HOME/miniconda3/etc/profile.d/conda.sh"

# ---------------------------------------------------------------------------
# Inputs and reusable source artifacts
# ---------------------------------------------------------------------------
REF_PDB="${REF_PDB:-$ROOT_DIR/raw/assets/structures/2B3P.pdb}"
REF_FASTA="${REF_FASTA:-$ROOT_DIR/raw/current_project/competition_package/AAseqs of 5 GFP proteins.txt}"
EXCLUSION_LIST="${EXCLUSION_LIST:-$ROOT_DIR/raw/current_project/competition_package/Exclusion_List.csv}"

RUN05_LOOP_CSV="${RUN05_LOOP_CSV:-$ROOT_DIR/run05/output/loop_plddt_analysis.csv}"
RUN05_PLDDT_MANIFEST="${RUN05_PLDDT_MANIFEST:-$ROOT_DIR/run05/top48_plddt_boltz_pdbs/manifest.tsv}"
RUN07_PSSM_DIR="${RUN07_PSSM_DIR:-$ROOT_DIR/run07/pssm}"
SOURCE_CONSERVATIVE_POSITIONS="${SOURCE_CONSERVATIVE_POSITIONS:-$RUN07_PSSM_DIR/conservative_positions.csv}"
SOURCE_MATRIX_CSV="${SOURCE_MATRIX_CSV:-$RUN07_PSSM_DIR/2b3p_sequence_variance_matrix.csv}"
SOURCE_POSITION_TIERS_CSV="${SOURCE_POSITION_TIERS_CSV:-$RUN07_PSSM_DIR/2b3p_position_tiers.csv}"
SOURCE_BIAS_JSON="${SOURCE_BIAS_JSON:-$RUN07_PSSM_DIR/ligandmpnn_bias_AA_per_residue_conservative.json}"
SOURCE_OMIT_JSON="${SOURCE_OMIT_JSON:-$RUN07_PSSM_DIR/ligandmpnn_omit_AA_per_residue_conservative.json}"

TRIMMED_BACKBONES_TSV="${TRIMMED_BACKBONES_TSV:-$INPUT_DIR/trimmed_backbones.tsv}"
INPUT_BACKBONES_TXT="${INPUT_BACKBONES_TXT:-$INPUT_DIR/input_backbones.txt}"
REDESIGNED_RESIDUES_TXT="${REDESIGNED_RESIDUES_TXT:-$PSSM_DIR/redesigned_residues_barrel_pssm.txt}"
BARREL_POSITIONS_CSV="${BARREL_POSITIONS_CSV:-$PSSM_DIR/barrel_pssm_positions.csv}"
BIAS_JSON="${BIAS_JSON:-$PSSM_DIR/ligandmpnn_bias_AA_per_residue_barrel_pssm.json}"
OMIT_JSON="${OMIT_JSON:-$PSSM_DIR/ligandmpnn_omit_AA_per_residue_barrel_pssm.json}"
MATRIX_CSV="${MATRIX_CSV:-$PSSM_DIR/2b3p_sequence_variance_matrix.csv}"

LIGANDMPNN_RUNPY="${LIGANDMPNN_RUNPY:-$HOME/Documents/packages/LigandMPNN/run.py}"
LIGANDMPNN_CKPT="${LIGANDMPNN_CKPT:-$HOME/Documents/packages/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt}"
TMALIGN_BIN="${TMALIGN_BIN:-$HOME/Documents/packages/TMalign}"

MPNN_OUT="${MPNN_OUT:-$OUT_DIR/ligandmpnn_barrel_pssm}"
RAW_FASTA="${RAW_FASTA:-$OUT_DIR/raw_mpnn_sequences.fasta}"
PROCESSED_FASTA="${PROCESSED_FASTA:-$OUT_DIR/processed_sequences.fasta}"
PREFILTER_CSV="${PREFILTER_CSV:-$OUT_DIR/sequence_prefilter.csv}"
PREFILTER_SUMMARY="${PREFILTER_SUMMARY:-$OUT_DIR/sequence_prefilter_summary.md}"
BOLTZ_FASTA="${BOLTZ_FASTA:-$OUT_DIR/boltz_candidates.fasta}"
BOLTZ_MANIFEST="${BOLTZ_MANIFEST:-$OUT_DIR/boltz_candidates_manifest.csv}"
BOLTZ_YAML_DIR="${BOLTZ_YAML_DIR:-$OUT_DIR/boltz_yamls}"
BOLTZ_CHUNK_DIR="${BOLTZ_CHUNK_DIR:-$OUT_DIR/boltz_yaml_chunks}"
BOLTZ_SMOKE_YAML_DIR="${BOLTZ_SMOKE_YAML_DIR:-$OUT_DIR/boltz_smoke_yamls}"
BOLTZ_SMOKE_OUT="${BOLTZ_SMOKE_OUT:-$OUT_DIR/boltz_smoke_results}"
BOLTZ_OUT="${BOLTZ_OUT:-$OUT_DIR/boltz_results}"
BOLTZ_ANALYSIS_CSV="${BOLTZ_ANALYSIS_CSV:-$OUT_DIR/boltz_analysis.csv}"
TOP_FASTA="${TOP_FASTA:-$OUT_DIR/top_candidates_by_plddt.fasta}"
TOP_FOLD_GATED_FASTA="${TOP_FOLD_GATED_FASTA:-$OUT_DIR/top_fold_gated_candidates.fasta}"
TOP_PDB_DIR="${TOP_PDB_DIR:-$RUN_DIR/top64_boltz_pdbs}"

# ---------------------------------------------------------------------------
# Environments, scale, gates, and toggles
# ---------------------------------------------------------------------------
ANALYSIS_ENV="${ANALYSIS_ENV:-base}"
LIGANDMPNN_ENV="${LIGANDMPNN_ENV:-ligandmpnn_env}"
BOLTZ_ENV="${BOLTZ_ENV:-boltz_env}"
MPNN_CUDA_VISIBLE_DEVICES="${MPNN_CUDA_VISIBLE_DEVICES:-1}"

RUN09_TOP_LOOP_N="${RUN09_TOP_LOOP_N:-2}"
FINAL_LENGTH="${FINAL_LENGTH:-230}"
EXPECTED_EDITABLE_COUNT="${EXPECTED_EDITABLE_COUNT:-52}"

MPNN_TEMPERATURES="${MPNN_TEMPERATURES:-0.05 0.10 0.15 0.20 0.30}"
MPNN_BATCH_SIZE="${MPNN_BATCH_SIZE:-100}"
MPNN_BATCHES_PER_BACKBONE_TEMPERATURE="${MPNN_BATCHES_PER_BACKBONE_TEMPERATURE:-100}"
MPNN_SEED="${MPNN_SEED:-42}"

BOLTZ_TARGET_TOTAL="${BOLTZ_TARGET_TOTAL:-2000}"
BOLTZ_PER_LANE="${BOLTZ_PER_LANE:-200}"
BOLTZ_CHUNK_SIZE="${BOLTZ_CHUNK_SIZE:-50}"
BOLTZ_SMOKE_N="${BOLTZ_SMOKE_N:-10}"
BOLTZ_BATCH_PAUSE_SECONDS="${BOLTZ_BATCH_PAUSE_SECONDS:-60}"
TOP_N_FOR_PACKAGE="${TOP_N_FOR_PACKAGE:-64}"

BOLTZ_CUDA_DEVICE_ORDER="${BOLTZ_CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
BOLTZ_5080_CUDA_VISIBLE_DEVICES="${BOLTZ_5080_CUDA_VISIBLE_DEVICES:-1}"
BOLTZ_5080_EXPECT_GPU_NAME="${BOLTZ_5080_EXPECT_GPU_NAME:-5080}"
BOLTZ_FALLBACK_CUDA_VISIBLE_DEVICES="${BOLTZ_FALLBACK_CUDA_VISIBLE_DEVICES:-0}"
BOLTZ_FALLBACK_EXPECT_GPU_NAME="${BOLTZ_FALLBACK_EXPECT_GPU_NAME:-A4000}"
BOLTZ_EXTRA_ARGS="${BOLTZ_EXTRA_ARGS:-}"

LM_CONF_MIN="${LM_CONF_MIN:-0.43}"
MUTATION_TARGET="${MUTATION_TARGET:-40}"
DIVERSITY_MAX_JACCARD="${DIVERSITY_MAX_JACCARD:-0.90}"
HOMOPOLYMER_MAX="${HOMOPOLYMER_MAX:-4}"
DIPEPTIDE_REPEAT_MAX="${DIPEPTIDE_REPEAT_MAX:-3}"
LOW_COMPLEXITY_MAX="${LOW_COMPLEXITY_MAX:-0.05}"
NET_CHARGE_MIN="${NET_CHARGE_MIN:--20}"
NET_CHARGE_MAX="${NET_CHARGE_MAX:-20}"
CYSTEINE_MAX="${CYSTEINE_MAX:-3}"

PLDDT_GATE="${PLDDT_GATE:-0.85}"
PLDDT_PREFERRED="${PLDDT_PREFERRED:-0.90}"
TM2_GATE="${TM2_GATE:-0.75}"
ALIGNED_GATE="${ALIGNED_GATE:-210}"
RMSD_GATE="${RMSD_GATE:-3.5}"

RUN09_DO_PREP="${RUN09_DO_PREP:-1}"
RUN09_DO_MPNN="${RUN09_DO_MPNN:-1}"
RUN09_DO_MERGE="${RUN09_DO_MERGE:-1}"
RUN09_DO_PREFILTER="${RUN09_DO_PREFILTER:-1}"
RUN09_DO_YAMLS="${RUN09_DO_YAMLS:-1}"
RUN09_DO_SMOKE="${RUN09_DO_SMOKE:-1}"
RUN09_DO_BOLTZ="${RUN09_DO_BOLTZ:-1}"
RUN09_DO_ANALYZE="${RUN09_DO_ANALYZE:-1}"
RUN09_DO_PACKAGE="${RUN09_DO_PACKAGE:-1}"
RUN09_CLEAN_MPNN_BULK="${RUN09_CLEAN_MPNN_BULK:-1}"
FORCE="${FORCE:-0}"

mkdir -p "$INPUT_DIR" "$PSSM_DIR" "$OUT_DIR" "$MPNN_OUT" "$BOLTZ_OUT"

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

temperature_tag() {
    local t="$1"
    echo "T${t//./p}"
}

cleanup_ligandmpnn_bulk() {
    local out_dir="$1"
    if [[ "$RUN09_CLEAN_MPNN_BULK" != "1" ]]; then
        return 0
    fi
    [[ -d "$out_dir" ]] || return 0
    find "$out_dir" -type d \( -name "backbones" -o -name "packed" -o -name "packed_pdbs" -o -name "pdbs" \) -prune -exec rm -rf {} +
    find "$out_dir" -type f \( -name "*.pdb" -o -name "*.cif" \) ! -path "*/seqs/*" -delete
}

preflight_boltz_device() {
    local expected="$1"
    export CUDA_DEVICE_ORDER="$BOLTZ_CUDA_DEVICE_ORDER"
    python - "$expected" <<'PY'
import os
import sys
import torch

expected = sys.argv[1]
if not torch.cuda.is_available():
    raise SystemExit("ERROR: CUDA is not available in the Boltz environment")
name = torch.cuda.get_device_name(0)
capability = torch.cuda.get_device_capability(0)
print(
    f"Boltz resolved CUDA device: {name} capability={capability} "
    f"CUDA_DEVICE_ORDER={os.environ.get('CUDA_DEVICE_ORDER', '')} "
    f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES', '')}",
    flush=True,
)
if expected and expected.lower() not in name.lower():
    raise SystemExit(f"ERROR: expected GPU name containing {expected!r}, resolved {name!r}")
PY
}

run_boltz_dir() {
    local yaml_dir="$1"
    local out_dir="$2"
    local log_file="$3"
    local visible="$4"
    local expected="$5"
    activate_env "$BOLTZ_ENV"
    export CUDA_DEVICE_ORDER="$BOLTZ_CUDA_DEVICE_ORDER"
    export CUDA_VISIBLE_DEVICES="$visible"
    if ! preflight_boltz_device "$expected"; then
        return 1
    fi
    set +e
    boltz predict "$yaml_dir" --out_dir "$out_dir" $BOLTZ_EXTRA_ARGS 2>&1 | tee -a "$log_file"
    local status="${PIPESTATUS[0]}"
    set -e
    return "$status"
}

chunk_complete() {
    local yaml_dir="$1"
    local out_dir="$2"
    local result_root="$out_dir/boltz_results_$(basename "$yaml_dir")"
    local yaml name pred_dir
    shopt -s nullglob
    for yaml in "$yaml_dir"/*.yaml; do
        name="$(basename "$yaml" .yaml)"
        pred_dir="$result_root/predictions/$name"
        if [[ ! -s "$pred_dir/${name}_model_0.cif" && ! -s "$pred_dir/${name}_model_0.pdb" ]]; then
            shopt -u nullglob
            return 1
        fi
    done
    shopt -u nullglob
    return 0
}

echo "Run09 loop-seeded barrel PSSM refinement workflow"
echo "Root:              $ROOT_DIR"
echo "Run dir:           $RUN_DIR"
echo "Trim sfGFP coord:  $FINAL_LENGTH"
echo "Top Run05 loops:   $RUN09_TOP_LOOP_N"
echo "Editable residues: expected $EXPECTED_EDITABLE_COUNT"
echo "Temperatures:      $MPNN_TEMPERATURES"
echo "MPNN scale/lane:   $MPNN_BATCH_SIZE x $MPNN_BATCHES_PER_BACKBONE_TEMPERATURE = $((MPNN_BATCH_SIZE * MPNN_BATCHES_PER_BACKBONE_TEMPERATURE))"
echo "Boltz target:      $BOLTZ_TARGET_TOTAL candidates; chunks of $BOLTZ_CHUNK_SIZE; pause ${BOLTZ_BATCH_PAUSE_SECONDS}s"
echo "Boltz 5080:        CUDA_DEVICE_ORDER=$BOLTZ_CUDA_DEVICE_ORDER CUDA_VISIBLE_DEVICES=$BOLTZ_5080_CUDA_VISIBLE_DEVICES expect=$BOLTZ_5080_EXPECT_GPU_NAME"
echo "Boltz fallback:    CUDA_VISIBLE_DEVICES=$BOLTZ_FALLBACK_CUDA_VISIBLE_DEVICES expect=$BOLTZ_FALLBACK_EXPECT_GPU_NAME"

# ---------------------------------------------------------------------------
# Step 1: Prepare trimmed Run05 loop backbones and barrel PSSM guidance
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_PREP" == "1" ]]; then
    log_step "Step 1: Prepare sfGFP-coordinate 1-${FINAL_LENGTH} loop-seeded backbones and barrel PSSM guidance"
    activate_env "$ANALYSIS_ENV"
    python "$SCRIPT_DIR/run09_loop_pssm_utils.py" prepare-inputs \
        --root-dir "$ROOT_DIR" \
        --loop-csv "$RUN05_LOOP_CSV" \
        --manifest "$RUN05_PLDDT_MANIFEST" \
        --ref-fasta "$REF_FASTA" \
        --conservative-positions "$SOURCE_CONSERVATIVE_POSITIONS" \
        --source-matrix "$SOURCE_MATRIX_CSV" \
        --source-position-tiers "$SOURCE_POSITION_TIERS_CSV" \
        --source-bias "$SOURCE_BIAS_JSON" \
        --source-omit "$SOURCE_OMIT_JSON" \
        --input-dir "$INPUT_DIR" \
        --pssm-dir "$PSSM_DIR" \
        --top-loop-n "$RUN09_TOP_LOOP_N" \
        --final-length "$FINAL_LENGTH" \
        --expected-editable-count "$EXPECTED_EDITABLE_COUNT"
else
    echo "Skipping prep (RUN09_DO_PREP=$RUN09_DO_PREP)"
fi

# ---------------------------------------------------------------------------
# Step 2: LigandMPNN sequence generation
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_MPNN" == "1" ]]; then
    log_step "Step 2: LigandMPNN beta-barrel PSSM refinement"
    activate_env "$LIGANDMPNN_ENV"
    export CUDA_VISIBLE_DEVICES="$MPNN_CUDA_VISIBLE_DEVICES"
    echo "LigandMPNN batch size: $MPNN_BATCH_SIZE"
    echo "LigandMPNN batches per backbone-temperature: $MPNN_BATCHES_PER_BACKBONE_TEMPERATURE"

    tail -n +2 "$TRIMMED_BACKBONES_TSV" | while IFS=$'\t' read -r \
        backbone_id loop_rank donor_name source_pdb trimmed_pdb length trimmed_tail_removed \
        trim_boundary_sf_position trim_boundary_donor_index missing_sf_positions_1_to_trim \
        inherited_source_differences index_to_sf loop_129_147 loop_129_147_seq_start \
        loop_129_147_seq_end loop_188_198 loop_188_198_seq_start loop_188_198_seq_end \
        loop_209_216 loop_209_216_seq_start loop_209_216_seq_end redesigned_residues_txt \
        bias_json omit_json sequence; do
        [[ -z "$backbone_id" ]] && continue
        REDESIGNED_RESIDUES="$(cat "$redesigned_residues_txt")"
        echo "Backbone $backbone_id length=$length trim_sf=$trim_boundary_sf_position trim_donor_index=$trim_boundary_donor_index"
        echo "Redesigned residues for $backbone_id: $REDESIGNED_RESIDUES"
        for temp in $MPNN_TEMPERATURES; do
            tag="$(temperature_tag "$temp")"
            bb_out="$MPNN_OUT/$backbone_id/$tag"
            fa_path="$bb_out/seqs/$(basename "$trimmed_pdb" .pdb).fa"
            if [[ "$FORCE" != "1" && -s "$fa_path" ]]; then
                echo "Skipping existing LigandMPNN output for $backbone_id $temp: $fa_path"
                continue
            fi
            mkdir -p "$bb_out"
            echo "LigandMPNN on $backbone_id temperature=$temp"
            mpnn_cmd=(
                python "$LIGANDMPNN_RUNPY"
                --seed "$MPNN_SEED"
                --model_type ligand_mpnn
                --pdb_path "$trimmed_pdb"
                --out_folder "$bb_out"
                --redesigned_residues "$REDESIGNED_RESIDUES"
                --bias_AA_per_residue "$bias_json"
                --omit_AA_per_residue "$omit_json"
                --ligand_mpnn_use_side_chain_context 1
                --temperature "$temp"
                --batch_size "$MPNN_BATCH_SIZE"
                --number_of_batches "$MPNN_BATCHES_PER_BACKBONE_TEMPERATURE"
                --save_backbones 0
                --pack_side_chains 0
                --checkpoint_ligand_mpnn "$LIGANDMPNN_CKPT"
            )
            "${mpnn_cmd[@]}" 2>&1 | tee -a "$OUT_DIR/ligandmpnn_barrel_pssm.log"
            cleanup_ligandmpnn_bulk "$bb_out"
        done
    done
else
    echo "Skipping LigandMPNN (RUN09_DO_MPNN=$RUN09_DO_MPNN)"
fi

# ---------------------------------------------------------------------------
# Step 3: Merge LigandMPNN FASTA outputs
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_MERGE" == "1" ]]; then
    log_step "Step 3: Merge LigandMPNN FASTA outputs"
    activate_env "$ANALYSIS_ENV"
    python "$SCRIPT_DIR/run09_loop_pssm_utils.py" merge-mpnn \
        --mpnn-out "$MPNN_OUT" \
        --out-fasta "$RAW_FASTA"
else
    echo "Skipping merge (RUN09_DO_MERGE=$RUN09_DO_MERGE)"
fi

# ---------------------------------------------------------------------------
# Step 4: Sequence filters and 2,000-candidate Boltz selection
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_PREFILTER" == "1" ]]; then
    log_step "Step 4: Filter sequences and select Boltz candidates"
    activate_env "$ANALYSIS_ENV"
    python "$SCRIPT_DIR/run09_loop_pssm_utils.py" prefilter \
        --raw-fasta "$RAW_FASTA" \
        --processed-fasta "$PROCESSED_FASTA" \
        --prefilter-csv "$PREFILTER_CSV" \
        --summary-md "$PREFILTER_SUMMARY" \
        --boltz-fasta "$BOLTZ_FASTA" \
        --boltz-manifest "$BOLTZ_MANIFEST" \
        --metadata-tsv "$TRIMMED_BACKBONES_TSV" \
        --barrel-positions "$BARREL_POSITIONS_CSV" \
        --matrix-csv "$MATRIX_CSV" \
        --ref-fasta "$REF_FASTA" \
        --exclusion-list "$EXCLUSION_LIST" \
        --final-length "$FINAL_LENGTH" \
        --lm-conf-min "$LM_CONF_MIN" \
        --homopolymer-max "$HOMOPOLYMER_MAX" \
        --dipeptide-repeat-max "$DIPEPTIDE_REPEAT_MAX" \
        --low-complexity-max "$LOW_COMPLEXITY_MAX" \
        --net-charge-min "$NET_CHARGE_MIN" \
        --net-charge-max "$NET_CHARGE_MAX" \
        --cysteine-max "$CYSTEINE_MAX" \
        --boltz-target-total "$BOLTZ_TARGET_TOTAL" \
        --boltz-per-lane "$BOLTZ_PER_LANE" \
        --mutation-target "$MUTATION_TARGET" \
        --diversity-max-jaccard "$DIVERSITY_MAX_JACCARD"
else
    echo "Skipping prefilter (RUN09_DO_PREFILTER=$RUN09_DO_PREFILTER)"
fi

# ---------------------------------------------------------------------------
# Step 5: Generate Boltz YAMLs and chunks
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_YAMLS" == "1" ]]; then
    log_step "Step 5: Generate Boltz YAMLs and chunks"
    activate_env "$ANALYSIS_ENV"
    python "$SCRIPT_DIR/run09_loop_pssm_utils.py" make-yamls \
        --fasta "$BOLTZ_FASTA" \
        --yaml-dir "$BOLTZ_YAML_DIR" \
        --chunk-dir "$BOLTZ_CHUNK_DIR" \
        --smoke-dir "$BOLTZ_SMOKE_YAML_DIR" \
        --chunk-size "$BOLTZ_CHUNK_SIZE" \
        --smoke-n "$BOLTZ_SMOKE_N"
else
    echo "Skipping YAML generation (RUN09_DO_YAMLS=$RUN09_DO_YAMLS)"
fi

# ---------------------------------------------------------------------------
# Step 6: Boltz-2 validation, 5080 first with A4000 fallback
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_BOLTZ" == "1" ]]; then
    log_step "Step 6: Boltz-2 prediction"
    fallback_needed=0
    if [[ "$RUN09_DO_SMOKE" == "1" ]]; then
        echo "Running 5080 smoke test on $BOLTZ_SMOKE_N YAMLs"
        rm -rf "$BOLTZ_SMOKE_OUT"
        if ! run_boltz_dir "$BOLTZ_SMOKE_YAML_DIR" "$BOLTZ_SMOKE_OUT" "$OUT_DIR/boltz_5080_smoke.log" "$BOLTZ_5080_CUDA_VISIBLE_DEVICES" "$BOLTZ_5080_EXPECT_GPU_NAME"; then
            echo "WARNING: 5080 smoke failed; switching full run to fallback GPU"
            fallback_needed=1
        fi
    fi

    if [[ "$fallback_needed" == "0" ]]; then
        first_chunk=1
        for chunk in "$BOLTZ_CHUNK_DIR"/chunk_*; do
            [[ -d "$chunk" ]] || continue
            if [[ "$first_chunk" == "0" && "$BOLTZ_BATCH_PAUSE_SECONDS" != "0" ]]; then
                echo "Pausing ${BOLTZ_BATCH_PAUSE_SECONDS}s before next 5080 Boltz chunk"
                sleep "$BOLTZ_BATCH_PAUSE_SECONDS"
            fi
            first_chunk=0
            if chunk_complete "$chunk" "$BOLTZ_OUT"; then
                echo "Skipping complete Boltz chunk on 5080: $(basename "$chunk")"
                continue
            fi
            echo "Running Boltz on 5080 for $(basename "$chunk")"
            if ! run_boltz_dir "$chunk" "$BOLTZ_OUT" "$OUT_DIR/boltz_5080_chunks.log" "$BOLTZ_5080_CUDA_VISIBLE_DEVICES" "$BOLTZ_5080_EXPECT_GPU_NAME"; then
                echo "WARNING: 5080 Boltz failed on $(basename "$chunk"); switching to A4000 fallback"
                fallback_needed=1
                break
            fi
        done
    fi

    if [[ "$fallback_needed" == "1" ]]; then
        first_chunk=1
        for chunk in "$BOLTZ_CHUNK_DIR"/chunk_*; do
            [[ -d "$chunk" ]] || continue
            if [[ "$first_chunk" == "0" && "$BOLTZ_BATCH_PAUSE_SECONDS" != "0" ]]; then
                echo "Pausing ${BOLTZ_BATCH_PAUSE_SECONDS}s before next fallback Boltz chunk"
                sleep "$BOLTZ_BATCH_PAUSE_SECONDS"
            fi
            first_chunk=0
            if chunk_complete "$chunk" "$BOLTZ_OUT"; then
                echo "Skipping complete Boltz chunk on fallback: $(basename "$chunk")"
                continue
            fi
            echo "Running/resuming Boltz on fallback GPU for $(basename "$chunk")"
            run_boltz_dir "$chunk" "$BOLTZ_OUT" "$OUT_DIR/boltz_a4000_fallback.log" "$BOLTZ_FALLBACK_CUDA_VISIBLE_DEVICES" "$BOLTZ_FALLBACK_EXPECT_GPU_NAME"
        done
    fi
else
    echo "Skipping Boltz prediction (RUN09_DO_BOLTZ=$RUN09_DO_BOLTZ)"
fi

# ---------------------------------------------------------------------------
# Step 7: Analyze Boltz predictions
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_ANALYZE" == "1" ]]; then
    log_step "Step 7: Analyze Boltz predictions"
    activate_env "$ANALYSIS_ENV"
    python "$SCRIPT_DIR/run09_loop_pssm_utils.py" analyze \
        --boltz-out "$BOLTZ_OUT" \
        --fasta "$BOLTZ_FASTA" \
        --candidate-manifest "$BOLTZ_MANIFEST" \
        --out-csv "$BOLTZ_ANALYSIS_CSV" \
        --top-fasta "$TOP_FASTA" \
        --fold-fasta "$TOP_FOLD_GATED_FASTA" \
        --ref-pdb "$REF_PDB" \
        --tmalign-bin "$TMALIGN_BIN" \
        --top-n "$TOP_N_FOR_PACKAGE" \
        --plddt-gate "$PLDDT_GATE" \
        --plddt-preferred "$PLDDT_PREFERRED" \
        --tm2-gate "$TM2_GATE" \
        --aligned-gate "$ALIGNED_GATE" \
        --rmsd-gate "$RMSD_GATE"
else
    echo "Skipping analysis (RUN09_DO_ANALYZE=$RUN09_DO_ANALYZE)"
fi

# ---------------------------------------------------------------------------
# Step 8: Package top structures
# ---------------------------------------------------------------------------
if [[ "$RUN09_DO_PACKAGE" == "1" ]]; then
    log_step "Step 8: Package top structures"
    activate_env "$ANALYSIS_ENV"
    python "$SCRIPT_DIR/run09_loop_pssm_utils.py" package \
        --analysis-csv "$BOLTZ_ANALYSIS_CSV" \
        --out-dir "$TOP_PDB_DIR" \
        --top-n "$TOP_N_FOR_PACKAGE"
else
    echo "Skipping package (RUN09_DO_PACKAGE=$RUN09_DO_PACKAGE)"
fi

log_step "Run09 workflow complete"
echo "Main outputs:"
echo "  Trimmed backbones: $TRIMMED_BACKBONES_TSV"
echo "  Barrel PSSM set:   $BARREL_POSITIONS_CSV"
echo "  Raw MPNN FASTA:    $RAW_FASTA"
echo "  Prefilter CSV:     $PREFILTER_CSV"
echo "  Boltz FASTA:       $BOLTZ_FASTA"
echo "  Boltz analysis:    $BOLTZ_ANALYSIS_CSV"
echo "  Top FASTA:         $TOP_FASTA"
echo "  Top PDB bundle:    $TOP_PDB_DIR"
