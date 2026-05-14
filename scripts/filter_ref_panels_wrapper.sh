#!/bin/bash
# filter_ref_panels_wrapper.sh
# Wrapper to load python and PLINK modules before filtering reference panels

set -e

SNP_LIST="$1"
REF_PANEL="$2"
CHROM="$3"
OUTPUT_PREFIX="$4"

if [[ ! -f "$SNP_LIST" ]]; then
    echo "Error: SNP list not found: $SNP_LIST" >&2
    exit 1
fi

module load slurm
module load python/3.10.2

# Load PLINK module from software config (fallback to known default)
PLINK_MODULE=$(python3 -c "
import sys
sys.path.insert(0, 'utils')
from bioconfigme import get_software_module
print(get_software_module('plink2'))
" 2>/dev/null)

if [[ -z "$PLINK_MODULE" ]]; then
    PLINK_MODULE="plink2/2.00a3.3lm"
fi

echo "Loading PLINK module: $PLINK_MODULE"
module load "$PLINK_MODULE"

python3 scripts/filter_ref_panels.py \
    --snp-list "$SNP_LIST" \
    --ref-panel "$REF_PANEL" \
    --chrom "$CHROM" \
    --output "$OUTPUT_PREFIX"
