# =============================================================================
# 01 | Global genetic correlation -- LDSC
# =============================================================================
# Estimates genome-wide genetic correlation (rg) between each MDD
# phenotype and each cardiovascular trait using cross-trait LD Score
# Regression (LDSC; Bulik-Sullivan et al., Nat Genet 2015).
#
# LDSC is an open-source Python command-line tool. This script is a fully
# transparent wrapper that (1) converts each GWAS to the .sumstats.gz format
# with munge_sumstats.py and (2) runs ldsc.py --rg for every MDD x CVD
# pair. No in-house package is required.
#
# One-time setup:
#   - Install LDSC:  https://github.com/bulik/ldsc
#   - Download the EUR reference files linked from the LDSC wiki:
#       reference/ldsc/w_hm3.snplist
#       reference/ldsc/eur_w_ld_chr/
# =============================================================================

source("config.R")

munge_dir <- file.path(RESULTS_DIR, "ldsc", "munged")
rg_dir    <- file.path(RESULTS_DIR, "ldsc", "rg")
dir.create(munge_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rg_dir,    recursive = TRUE, showWarnings = FALSE)

all_traits <- c(MDD_TRAITS, CVD_TRAITS)

## ---- Step 1: munge each GWAS into the LDSC .sumstats format -----------------
for (trait in all_traits) {
  out <- file.path(munge_dir, trait)
  if (file.exists(paste0(out, ".sumstats.gz"))) next   # skip if already done
  system2(LDSC_PY, c(
    file.path(LDSC_DIR, "munge_sumstats.py"),
    "--sumstats",       file.path(DATA_DIR, paste0(trait, ".txt")),
    "--out",            out,
    "--merge-alleles",  file.path(REF_DIR, "ldsc", "w_hm3.snplist"),
    "--snp",            COL$SNP,
    "--a1",             COL$EA,
    "--a2",             COL$OA,
    "--p",              COL$P,
    "--N-col",          COL$N,
    "--frq",            COL$EAF,
    "--signed-sumstats", paste0(COL$BETA, ",0"),
    "--chunksize",      "500000"
  ))
}

## ---- Step 2: cross-trait rg for every MDD x CVD pair -------------------
ld_chr <- file.path(REF_DIR, "ldsc", "eur_w_ld_chr/")

for (m in MDD_TRAITS) {
  for (cv in CVD_TRAITS) {
    system2(LDSC_PY, c(
      file.path(LDSC_DIR, "ldsc.py"),
      "--rg", paste0(file.path(munge_dir, m),  ".sumstats.gz,",
                     file.path(munge_dir, cv), ".sumstats.gz"),
      "--ref-ld-chr", ld_chr,
      "--w-ld-chr",   ld_chr,
      "--out", file.path(rg_dir, paste0(m, "_", cv, "_ldsc_rg"))
    ))
  }
}

cat("LDSC genetic correlation complete. Results in:", rg_dir, "\n")
cat("rg, SE and P are reported in the 'Summary of Genetic Correlation",
    "Results' block of each .log file.\n")