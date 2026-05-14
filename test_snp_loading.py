#!/usr/bin/env python3
"""Quick test to verify SNP loading from different GWAS file formats."""

import sys
import os

# Add scripts to path to import functions
sys.path.insert(0, "scripts")

# Import just the SNP loading function
import importlib.util
spec = importlib.util.spec_from_file_location("summarize", "scripts/summarize_dataset_loci_with_genes.py")
module = importlib.util.module_from_spec(spec)
sys.modules["summarize"] = module
spec.loader.exec_module(module)

# Test SNP loading
print("Testing SNP loading from input/ directory...")
loci_snps = module.load_snps_for_loci()

print(f"\n=== SNP Loading Results ===")
print(f"Total datasets loaded: {len(loci_snps)}")
for dataset, snps_df in sorted(loci_snps.items()):
    if snps_df is not None:
        print(f"  {dataset}: {len(snps_df)} SNPs")
        print(f"    CHR range: {snps_df['CHR'].unique()}")
        print(f"    POS range: {snps_df['POS'].min()}-{snps_df['POS'].max()}")

if len(loci_snps) == 0:
    print("[ERROR] No SNPs loaded!")
    sys.exit(1)
else:
    print("\n[SUCCESS] SNP loading working!")
    sys.exit(0)
