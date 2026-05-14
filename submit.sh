#!/bin/bash
#SBATCH --job-name=gwas_finemap
#SBATCH --output=logs/gwas_finemap_%j.out
#SBATCH --error=logs/gwas_finemap_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8

# Create logs directory
mkdir -p logs

# Parse command line arguments
DRY_RUN=""
TARGETS=()
EXTRA_FLAGS=""
SNAKEFILE=""
CORES=""
JOBS=""
UNTIL_RULE=""
ALLOWED_RULES=()

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run|-n)
            DRY_RUN="--dry-run"
            shift
            ;;
        --until)
            UNTIL_RULE="--until $2"
            shift 2
            ;;
        --allowed-rules)
            IFS=',' read -r -a _rules <<< "$2"
            for _r in "${_rules[@]}"; do
                if [[ -n "$_r" ]]; then
                    ALLOWED_RULES+=("$_r")
                fi
            done
            shift 2
            ;;
        --snakefile)
            SNAKEFILE="$2"
            shift 2
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --jobs|-j)
            JOBS="$2"
            shift 2
            ;;
        --*)
            # Pass through other flags
            EXTRA_FLAGS="$EXTRA_FLAGS $1"
            if [[ $2 && ! $2 =~ ^-- ]]; then
                EXTRA_FLAGS="$EXTRA_FLAGS $2"
                shift 2
            else
                shift
            fi
            ;;
        *)
            # This is a target file
            TARGETS+=("$1")
            shift
            ;;
    esac
done

# Detect if this is an individual rule submission
INDIVIDUAL_RULE=false
if [[ -n "$SNAKEFILE" && "$SNAKEFILE" =~ ^rules/ ]]; then
    INDIVIDUAL_RULE=true
fi

# Set defaults based on submission type
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    # Individual rule defaults
    DEFAULT_SNAKEFILE="$SNAKEFILE"
    DEFAULT_CORES="${CORES:-2}"
    DEFAULT_JOBS=1
    DEFAULT_MEM="128G"
    DEFAULT_TIME="2:00:00"
    DEFAULT_CPUS=2
else
    # Full pipeline defaults
    DEFAULT_SNAKEFILE="${SNAKEFILE:-Snakefile}"
    DEFAULT_CORES="${CORES:-8}"
    DEFAULT_JOBS=15
    DEFAULT_MEM="128G"
    DEFAULT_TIME="2-00:00:00"
    DEFAULT_CPUS=8
fi

# Default to 'all' if no targets specified and not individual rule
if [ ${#TARGETS[@]} -eq 0 ] && [[ "$INDIVIDUAL_RULE" == false ]]; then
    TARGETS=("all")
fi

# Set up environment
module load slurm
module load python/3.10.2

# Export Python path if needed
export PYTHONPATH="/s/nath-lab/alsamman/____MyCodes____/FineMappingSuite:$PYTHONPATH"

# Ensure snakemake state dirs exist before execution (helps on shared/NFS filesystems)
mkdir -p .snakemake .snakemake/log .snakemake/locks .snakemake/metadata .snakemake/incomplete

# Older Snakemake versions can race on metadata bookkeeping for clustered jobs.
# Use --drop-metadata when available to avoid false failures after successful jobs.
METADATA_FLAG=""
if snakemake --help 2>/dev/null | grep -q -- "--drop-metadata"; then
    METADATA_FLAG="--drop-metadata"
fi

# Print submission info
echo "Submission type: $([ "$INDIVIDUAL_RULE" == true ] && echo "Individual Rule" || echo "Full Pipeline")"
echo "Snakefile: $DEFAULT_SNAKEFILE"
echo "Cores: $DEFAULT_CORES"
echo "Target(s): ${TARGETS[*]:-"(from snakefile)"}"
if [[ ${#ALLOWED_RULES[@]} -gt 0 ]]; then
    echo "Allowed rule(s): ${ALLOWED_RULES[*]}"
fi
echo "Jobs: $DEFAULT_JOBS"
echo ""

ALLOWED_RULES_ARGS=()
if [[ ${#ALLOWED_RULES[@]} -gt 0 ]]; then
    ALLOWED_RULES_ARGS=(--allowed-rules "${ALLOWED_RULES[@]}")
fi

TARGET_ARGS=()
if [[ ${#TARGETS[@]} -gt 0 ]]; then
    TARGET_ARGS=(-- "${TARGETS[@]}")
fi

# Adjust SLURM header for individual rules
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    echo "Adjusting SLURM parameters for individual rule execution..."
fi

# Submit the workflow to SLURM
if [[ "$INDIVIDUAL_RULE" == true ]]; then
    # Individual rule submission - simpler configuration
    snakemake \
        --cluster "sbatch \
            --mem=$DEFAULT_MEM \
            --time=$DEFAULT_TIME \
            --job-name=rule_{rule}_{wildcards} \
            --cpus-per-task=$DEFAULT_CPUS \
            --output=logs/{rule}_{wildcards}_%j.out \
            --error=logs/{rule}_{wildcards}_%j.err" \
        --jobs ${JOBS:-$DEFAULT_JOBS} \
        --latency-wait 60 \
        --keep-going \
        --rerun-incomplete \
        --configfile configs/analysis.yml \
        --snakefile "$DEFAULT_SNAKEFILE" \
        --cores $DEFAULT_CORES \
        --verbose \
        --printshellcmds \
        --stats logs/snakemake_individual_stats.json \
        $METADATA_FLAG \
        $DRY_RUN \
        $UNTIL_RULE \
        "${ALLOWED_RULES_ARGS[@]}" \
        $EXTRA_FLAGS \
        "${TARGET_ARGS[@]}"
    SNAKEMAKE_EXIT_CODE=$?
else
    # Full pipeline submission - original configuration
    snakemake \
        --cluster "sbatch \
            --mem=128G \
            --time={resources.time} \
            --job-name={rule}_{wildcards} \
            --cpus-per-task={threads} \
            --output=logs/{rule}_{wildcards}_%j.out \
            --error=logs/{rule}_{wildcards}_%j.err" \
        --jobs $DEFAULT_JOBS \
        --latency-wait 120 \
        --keep-going \
        --rerun-incomplete \
        --configfile configs/analysis.yml \
        --snakefile "$DEFAULT_SNAKEFILE" \
        --cores $DEFAULT_CORES \
        --verbose \
        --printshellcmds \
        --stats logs/snakemake_stats.json \
        $METADATA_FLAG \
        $DRY_RUN \
        $UNTIL_RULE \
        "${ALLOWED_RULES_ARGS[@]}" \
        $EXTRA_FLAGS \
        "${TARGET_ARGS[@]}"
    SNAKEMAKE_EXIT_CODE=$?
fi

# Print completion message
if [[ ${SNAKEMAKE_EXIT_CODE:-1} -eq 0 ]]; then
    echo "Snakemake workflow completed for target(s): ${TARGETS[*]}"
    echo "Monitor jobs with: squeue -u $USER"
else
    echo "Snakemake workflow failed for target(s): ${TARGETS[*]}" >&2
fi

exit ${SNAKEMAKE_EXIT_CODE:-1}