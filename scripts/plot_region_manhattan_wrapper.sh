#!/bin/bash
# plot_region_manhattan_wrapper.sh
# Wrapper script to call the region Manhattan plotting R script using the configured R module.

set -e

INPUT_FILE="$1"
OUTPUT_PNG="$2"
REGION_NAME="$3"
REGION_CHR="$4"
REGION_START="$5"
REGION_END="$6"
DATASET="$7"
TABLE="$8"
TABLE_COLUMNS="$9"
SCALING_FILE="${10}"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file not found: $INPUT_FILE" >&2
    exit 1
fi

if [[ -z "$SCALING_FILE" || ! -f "$SCALING_FILE" ]]; then
    echo "Error: Global scaling file not found: $SCALING_FILE" >&2
    exit 1
fi

echo "Creating region Manhattan plot for $REGION_NAME ($DATASET - $TABLE)"
echo "Input: $INPUT_FILE"
echo "Output PNG: $OUTPUT_PNG"
echo "Scaling: $SCALING_FILE"

# Load required modules in the current shell
module load slurm
module load python/3.10.2

# Load R module from software.yml using bioconfigme
R_MODULE=$(python3 -c "
import sys
sys.path.insert(0, 'utils')
from bioconfigme import get_software_module
print(get_software_module('r'))
" 2>/dev/null)

if [[ -z "$R_MODULE" ]]; then
    echo "Warning: Could not load R module from bioconfigme, using default"
    R_MODULE="R/4.4.1-mkl"
fi

echo "Loading R module: $R_MODULE"
module load "$R_MODULE" 2>/dev/null || echo "Warning: R module load failed, trying direct path"

# Execute R script
Rscript scripts/plot_region_manhattan.R \
    "$INPUT_FILE" \
    "$OUTPUT_PNG" \
    "$REGION_NAME" \
    "$REGION_CHR" \
    "$REGION_START" \
    "$REGION_END" \
    "$DATASET" \
    "$TABLE" \
    "$TABLE_COLUMNS" \
    "$SCALING_FILE"

if [[ $? -eq 0 ]]; then
    echo "Region Manhattan plot successfully created."
else
    echo "Error: Failed to create region Manhattan plot." >&2
    exit 1
fi