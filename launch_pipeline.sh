#!/bin/bash
#SBATCH --job-name=Mito_Pipeline_Manager         # Job name in the queue
#SBATCH --cpus-per-task=2                        # Request 2 CPU cores
#SBATCH --mem=8G                                 # Request 8 GB of memory
#SBATCH --time=24:00:00                          # Max runtime of 24 hours
#SBATCH --output=log/mito_pipeline_manager_%j.log # Log file name (%j = job ID)

# Exit immediately if any command fails
set -e

# ==============================================================================
#                         User Configuration (edit below)
# ==============================================================================

# 1) File and Directory Paths
SCRIPT_BASE_DIR="/path/to/human_mito_calling-main"   # Pipeline root
MAIN_NF_SCRIPT="$SCRIPT_BASE_DIR/main.nf"                                      # Main Nextflow script
NEXTFLOW_CONFIG="$SCRIPT_BASE_DIR/nextflow_Bouchet.config"                     # Nextflow config

MASTER_SAMPLE_LIST="/path/to/test_sample_cram.tsv"  # TSV with all samples/URLs
OUTPUT_DIR="/path/to/output_dir"             # Final output dir (fallback)
WORK_DIR_BASE="/path/to/work_dir"                   # Base work dir for intermediates

# 2) Slurm Job Array Concurrency
CONCURRENT_SAMPLES=3  # Max number of array tasks to run simultaneously

# 3) Pipeline Mode: "population" or "disease"
PIPELINE_MODE="disease"

# Required only if PIPELINE_MODE="disease"; ignored for "population"
DISEASE_META_FILE="/path/to/disease_meta.tsv"

# 4) Environment for Nextflow (please change this per your server)
module load Nextflow/24.10.2

# ==============================================================================
#                                  Script Modes
# ------------------------------------------------------------------------------
# 1) Master Mode (default, no args, login node):
#       - Build unique sample list
#       - Submit Stage 1 (array workers) and Stage 2 (finalizer) jobs
# 2) Worker Mode (compute node; when $SLURM_ARRAY_TASK_ID is set):
#       - Process one unique sample (one array task)
# 3) Finalizer Mode (compute node; when $1 == --finalize):
#       - Merge/annotate final outputs across samples (Hail/Nextflow)
# ==============================================================================

# Internal file (do not change)
INTERNAL_SAMPLE_LIST="unique_samples_for_job_array.list"

# ==============================================================================
#                            Mode 3: Finalizer Mode
# ==============================================================================
if [ "$1" == "--finalize" ]; then
    echo "========================================================"
    echo "--- STAGE 2: FINALIZER MODE - Starting Nextflow Finalizer"
    echo "========================================================"

    MERGED_INPUT_DIR="${OUTPUT_DIR}/merged_results/wdl_output"

    if [ ! -d "$MERGED_INPUT_DIR" ]; then
        echo "[SKIP] 'merged_results' directory not found: $MERGED_INPUT_DIR"
        echo "[SKIP] No passing samples; skipping Stage 2."
        exit 0
    fi

    has_files=$(find "$MERGED_INPUT_DIR" -maxdepth 1 -type f \
        \( -name "*.merged.final.split.vcf" \
           -o -name "*.merged.per_base_coverage.tsv" \
           -o -name "*.merged.haplocheck_contamination.txt" \) -print -quit || true)

    if [ -z "$has_files" ]; then
        echo "[SKIP] No final artifacts detected under: $MERGED_INPUT_DIR"
        echo "[SKIP] Skipping Stage 2 (no annotation outputs to merge)."
        exit 0
    fi

    PROJECT_ROOT=$(dirname "$MAIN_NF_SCRIPT")
    FINALIZER_WORK_DIR="${WORK_DIR_BASE}/final"
    MERGE_NF_SCRIPT="$PROJECT_ROOT/merge.nf"

    if [ ! -f "$MERGE_NF_SCRIPT" ]; then
        echo "[ERROR] merge.nf not found at: $MERGE_NF_SCRIPT"
        exit 1
    fi

    echo "[*] Launching Nextflow finalizer: $MERGE_NF_SCRIPT"
    echo "[*] Input directory: $MERGED_INPUT_DIR"

    # Extra args for disease mode
    EXTRA_ARGS=()
    if [ "${PIPELINE_MODE}" = "disease" ]; then
        if [ -z "${DISEASE_META_FILE:-}" ] || [ ! -f "${DISEASE_META_FILE}" ]; then
            echo "[ERROR] PIPELINE_MODE='disease' but DISEASE_META_FILE is missing or not a file."
            echo "        Set DISEASE_META_FILE=/abs/path/to/disease_meta.tsv"
            exit 1
        fi
        echo "[*] Disease mode: using meta file ${DISEASE_META_FILE}"
        EXTRA_ARGS+=( --disease_meta_file "${DISEASE_META_FILE}" )
    fi

    nextflow run "$MERGE_NF_SCRIPT" \
        -c "$NEXTFLOW_CONFIG" \
        -profile cluster \
        -resume \
        --merged_dir "$MERGED_INPUT_DIR" \
        --pipeline_mode "$PIPELINE_MODE" \
        --outdir "$OUTPUT_DIR" \
        "${EXTRA_ARGS[@]}" \
        -w "$FINALIZER_WORK_DIR" \
        -with-report   "$WORK_DIR/running_report_final.html"

    NF_EXIT=$?
    if [ $NF_EXIT -eq 0 ]; then
        echo "[SUCCESS] Stage 2 (Nextflow finalizer) completed successfully."
    else
        echo "[ERROR] Stage 2 (Nextflow finalizer) failed with exit code $NF_EXIT."
        exit $NF_EXIT
    fi
    exit 0

# ==============================================================================
#                            Mode 2: Worker Mode
# ==============================================================================
elif [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    echo "================================================================="
    echo "--- STAGE 1: WORKER MODE on Compute Node (Task ${SLURM_ARRAY_TASK_ID}) ---"
    echo "================================================================="

    if [ ! -f "$INTERNAL_SAMPLE_LIST" ]; then
        echo "[ERROR] Unique sample list not found: $INTERNAL_SAMPLE_LIST"
        echo "[ERROR] This file should be created by Master Mode."
        exit 1
    fi

    # Map array index to the corresponding SAMPLE_ID (line = index + 1)
    TASK_INDEX=$((SLURM_ARRAY_TASK_ID + 1))
    SAMPLE_ID=$(awk -v line=$TASK_INDEX 'NR==line {print; exit}' "$INTERNAL_SAMPLE_LIST" | tr -d '\r')

    if [ -z "$SAMPLE_ID" ]; then
        echo "[ERROR] No SAMPLE_ID found in $INTERNAL_SAMPLE_LIST for task ${SLURM_ARRAY_TASK_ID} (line $TASK_INDEX)"
        exit 1
    fi
    echo "[*] Processing unique sample: ${SAMPLE_ID}"

    WORK_DIR="${WORK_DIR_BASE}/batch_${SLURM_ARRAY_TASK_ID}"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"
    echo "[*] Working directory: $(pwd)"

    # Build per-task input TSV containing all rows for this SAMPLE_ID
    TEMP_TSV="${WORK_DIR}/sample_${SAMPLE_ID}.tsv"
    awk -v id="$SAMPLE_ID" '$1==id {print}' "$MASTER_SAMPLE_LIST" > "$TEMP_TSV"

    echo "[*] Generated input file: ${TEMP_TSV}"
    echo "[*] Contents:"
    cat "${TEMP_TSV}"

    # Assemble Nextflow CLI arguments
    NEXTFLOW_RUN_ARGS=()
    NEXTFLOW_RUN_ARGS+=("-c" "$NEXTFLOW_CONFIG")
    NEXTFLOW_RUN_ARGS+=("-profile" "cluster")
    NEXTFLOW_RUN_ARGS+=("-resume")
    NEXTFLOW_RUN_ARGS+=("-w" "$WORK_DIR")
    NEXTFLOW_RUN_ARGS+=("--input" "$TEMP_TSV")
    NEXTFLOW_RUN_ARGS+=("--outdir" "$OUTPUT_DIR")
    NEXTFLOW_RUN_ARGS+=("-with-report" "$WORK_DIR/running_report_main.html")

    if [ "$PIPELINE_MODE" == "disease" ]; then
        echo "[*] Mode: disease (will check disease meta file)"
        if [ -z "$DISEASE_META_FILE" ] || [ ! -f "$DISEASE_META_FILE" ]; then
            echo "[ERROR] PIPELINE_MODE='disease' but DISEASE_META_FILE not set or missing."
            echo "[ERROR] Looked for: $DISEASE_META_FILE"
            exit 1
        fi
        NEXTFLOW_RUN_ARGS+=("--disease_meta_file" "$DISEASE_META_FILE")
    elif [ "$PIPELINE_MODE" == "population" ]; then
        echo "[*] Mode: population (disease meta not used)"
    else
        echo "[ERROR] Invalid PIPELINE_MODE: '$PIPELINE_MODE' (must be 'population' or 'disease')"
        exit 1
    fi

    echo "[*] Starting Nextflow..."
    nextflow run "$MAIN_NF_SCRIPT" "${NEXTFLOW_RUN_ARGS[@]}"

    NF_EXIT=$?
    if [ $NF_EXIT -eq 0 ]; then
        echo "Batch ${SLURM_ARRAY_TASK_ID} completed successfully for sample ${SAMPLE_ID}."
    else
        echo "Batch ${SLURM_ARRAY_TASK_ID} failed (exit code $NF_EXIT) for sample ${SAMPLE_ID}."
        echo "Retaining work directory for debugging: ${WORK_DIR}"
    fi
    echo "--- Finished Job Array Task ${SLURM_ARRAY_TASK_ID} ---"

# ==============================================================================
#                            Mode 1: Master Mode
# ==============================================================================
else
    echo "========================================================"
    echo "--- STAGE 0: MASTER MODE on Login Node - Submitting jobs"
    echo "========================================================"

    LOG_DIR="log"
    mkdir -p "$LOG_DIR"
    echo "[*] Log directory: $(pwd)/${LOG_DIR}"

    echo "[*] Building unique sample list from: $MASTER_SAMPLE_LIST"
    cut -f1 "$MASTER_SAMPLE_LIST" | grep -vi "^sample" | sort -u > "$INTERNAL_SAMPLE_LIST"

    NUM_SAMPLES=$(wc -l < "$INTERNAL_SAMPLE_LIST")
    if [ "$NUM_SAMPLES" -eq 0 ]; then
        echo "[ERROR] No valid sample IDs (after excluding header)."
        exit 1
    fi

    ARRAY_INDEX=$((NUM_SAMPLES - 1))
    echo "[*] Found ${NUM_SAMPLES} unique samples (saved to $INTERNAL_SAMPLE_LIST)."
    echo "[*] Submitting job array (0-${ARRAY_INDEX}) with concurrency ${CONCURRENT_SAMPLES}..."

    STAGE1_JOB_ID=$(sbatch --parsable --array=0-${ARRAY_INDEX}%${CONCURRENT_SAMPLES} \
        --output="${LOG_DIR}/mito_pipeline_worker_%A_%a.log" "$0")
    echo "[*] Stage 1 submitted (Job Array ID: ${STAGE1_JOB_ID})"

    STAGE2_JOB_ID=$(sbatch --parsable --dependency=afterok:${STAGE1_JOB_ID} \
        --output="${LOG_DIR}/mito_pipeline_finalizer_%A.log" "$0" --finalize)
    echo "[*] Stage 2 submitted (Job ID: ${STAGE2_JOB_ID}); will start after array ${STAGE1_JOB_ID} completes successfully."

    echo "[SUCCESS] All jobs submitted. Monitor with 'squeue -u $USER' or check the '${LOG_DIR}' directory."
fi
