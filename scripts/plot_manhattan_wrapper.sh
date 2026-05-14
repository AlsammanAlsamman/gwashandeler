#!/bin/bash
# plot_manhattan_wrapper.sh
# Wrapper script to call the R Manhattan plotting script

set -e

INPUT_FILE="$1"
OUTPUT_PNG="$2"
DATASET="$3"
TABLE="$4"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file not found: $INPUT_FILE" >&2
    exit 1
fi

echo "Creating Manhattan plot for $DATASET - $TABLE"
echo "Input: $INPUT_FILE"
echo "Output PNG: $OUTPUT_PNG"

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
Rscript scripts/plot_manhattan.R "$INPUT_FILE" "$OUTPUT_PNG" "$DATASET" "$TABLE"

if [[ $? -eq 0 ]]; then
    echo "Manhattan plot successfully created."
else
    echo "Error: Failed to create Manhattan plot." >&2
    exit 1
fi
