# Step 7 Loci Merging - REVIEW & FIXES APPLIED

## Rule Status: ✅ CORRECT

The Snakemake rule is **correctly configured** to:

### Input (Per-Dataset Loci from Step 6)
```
input: 
  - results/loci/AFR_AMR/AFR_AMR_loci.tsv
  - results/loci/AFR_EAS/AFR_EAS_loci.tsv
  - results/loci/AFR_EUR/AFR_EUR_loci.tsv
  - results/loci/AFR_gwas/AFR_gwas_loci.tsv
  - results/loci/AMR_EAS/AMR_EAS_loci.tsv
  - results/loci/AMR_EUR/AMR_EUR_loci.tsv
  - results/loci/AMR_gwas/AMR_gwas_loci.tsv
  - results/loci/BJEUR/BJEUR_loci.tsv
  - results/loci/EAS_gwas/EAS_gwas_loci.tsv
  - results/loci/EUR_EAS/EUR_EAS_loci.tsv
  - results/loci/EUR_gwas/EUR_gwas_loci.tsv
  - results/loci/Hisp/Hisp_loci.tsv
  - results/loci/WangEAS/WangEAS_loci.tsv
```
✅ Each is a **per-dataset merged loci TSV** from Step 6

### Processing (22 Parallel Jobs + Combine)
```
Rule: summarize_dataset_loci_with_genes_by_chr
  - Processes EACH CHROMOSOME (1-22) in parallel
  - Each job receives ALL 13 per-dataset loci TSVs
  - Filters to chromosome, merges with LD policy, outputs chr{N}.summary.tsv
  - 22 jobs run simultaneously (--jobs 22)

Rule: summarize_dataset_loci_with_genes
  - Collects all 22 chr TSVs
  - Merges, sorts by chr+pos
  - Generates final dataset_loci_gene_summary.tsv/.xlsx
```

## Critical Fixes Applied

### Fix 1: Multi-Format SNP Column Handling ✅
**Problem:** GWAS files have different column names
```
AFR/AMR/EAS/EUR.gwas.tsv:  CHR, POS, ID
WangEAS.tsv:              chrom, pos, varid / snpid / markername
Hisp.tsv:                 chrom, pos, rsid
```

**Solution:** Updated `load_snps_for_loci()` to:
- Try multiple column name variants for each field:
  - CHR variants: [CHR, CHROM, chrom, Chr, chromosome]
  - POS variants: [POS, POS_BUILD37, pos, Pos, position]
  - ID variants: [ID, SNP, SNPID, snpid, rsid, varid, markername]
- Auto-detect and map to standard format (CHR, POS, ID)
- Load max 1M SNPs per file (memory efficient)
- Validate all required columns present before processing

### Fix 2: PLINK Executable Discovery ✅
**Problem:** `[warn] PLINK executable not found in PATH`

**Solution:** Added `find_plink_executable()` function:
```python
def find_plink_executable():
    # Try PATH first
    plink_path = shutil.which("plink")
    if plink_path:
        return plink_path
    
    # Try common locations
    common_paths = [
        "/usr/bin/plink",
        "/usr/local/bin/plink",
        "/opt/plink/plink",
        os.path.expanduser("~/bin/plink"),
        os.path.expanduser("~/plink/plink"),
    ]
    for p in common_paths:
        if os.path.exists(p) and os.access(p, os.X_OK):
            return p
    
    return None
```

### Fix 3: Updated PLINK Call ✅
Now uses discovered executable instead of hardcoded "plink":
```python
cmd = [
    plink_exe,  # ← Uses found executable
    "--bfile", ref_panel_prefix,
    "--extract", all_file,
    "--r2",
    "--ld-window-kb", "10000",
    "--ld-window", "99999",
    "--ld-window-r2", "0",
    "--allow-no-sex",
    "--out", out_prefix,
]
```

## Data Flow Clarification

**NOT "merged loci" in the bad sense:**
```
Step 6 OUTPUT (per-dataset):
  └─ {dataset}_loci.tsv
     Contains: all loci for that dataset
     (already merged within that dataset's tables)

Step 7 TASK:
  ├─ Load loci from ALL 13 datasets
  ├─ Pool them together
  ├─ Merge intervals across ALL datasets
  │  (gap ≤ 100kb → always merge)
  │  (gap 150-250kb → LD test, merge if R² ≥ 0.2)
  │  (gap > 250kb → don't merge)
  ├─ Assign source loci to merged loci
  ├─ Annotate merged loci with genes
  ├─ Build per-dataset x table summary
  └─ Output final summary TSV/XLSX
```

This is **correct behavior** - Step 7 merges loci **across datasets** based on LD.

## Remaining TODO for PLINK

To get LD merging working:
1. **Verify PLINK is installed:**
   ```bash
   which plink
   # OR check:
   ls -la /usr/bin/plink /usr/local/bin/plink /opt/plink/plink
   ```

2. **If not found, install:**
   ```bash
   # Linux (Ubuntu/Debian)
   apt-get install plink
   
   # Or download from https://www.cog-genomics.org/plink/
   # Then add to PATH or one of the searched locations
   ```

3. **Verify PLINK works:**
   ```bash
   plink --version
   ```

## Configuration for LD Merging (analysis.yml)

```yaml
loci_identiffication:
  use_ld: true                          # Enable LD-based merging
  auto_merge_gap_kb: 100                # Auto-merge if gap ≤ 100kb
  ld_test_gap_min_kb: 150               # Run LD test if gap 150-250kb
  ld_test_gap_max_kb: 250
  ld_threshold: 0.2                     # Merge if R² ≥ 0.2
  n_boundary_snps: 10                   # Extract 10 SNPs from each side
  min_locus_width_kbp: 1                # Keep loci ≥ 1kb wide

ref_panels:
  EAS: /s/nath-lab/alsamman/____MyCodes____/resources/plink_by_superpop/1kg_phase3_hg19.EAS
  EUR: /s/nath-lab/alsamman/____MyCodes____/resources/plink_by_superpop/1kg_phase3_hg19.EUR
  AFR: /s/nath-lab/alsamman/____MyCodes____/resources/plink_by_superpop/1kg_phase3_hg19.AFR
  AMR: /s/nath-lab/alsamman/____MyCodes____/resources/plink_by_superpop/1kg_phase3_hg19.AMR
```

## File Changes Summary

✅ `scripts/summarize_dataset_loci_with_genes.py`:
  - Added `shutil` import
  - Added `find_plink_executable()` function
  - Updated `load_snps_for_loci()` for multi-format column detection
  - Updated `compute_ld_with_plink()` to use found executable
  - Fixed SNP loading debug output

✅ `rules/07_loci_dataset_summary.smk`:
  - Added `params.chromosome` to pass wildcard value
  - Added `wildcard_constraints: chromosome="\d+"`
  - Log directory creation in shell

✅ `steps.sh`:
  - No changes needed (already correct)
