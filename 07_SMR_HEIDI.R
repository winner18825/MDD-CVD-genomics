# =============================================================================
# 07 | SMR + HEIDI (multi-SNP)
# =============================================================================
# Summary-data-based Mendelian Randomization (Zhu et al., Nat Genet 2016)
# integrating GWAS summary statistics with GTEx v8 cis-eQTL data to identify
# genes whose expression is pleiotropically associated with each trait. The
# HEIDI test distinguishes genuine pleiotropy from linkage.
#
# SMR is an open-source command-line binary -- no in-house wrapper required:
#   https://yanglab.westlake.edu.cn/software/smr/
#
# Inputs (place under reference/ and data/):
#   - data/gwas/ma/<TRAIT>.ma          GWAS in SMR .ma format
#                                      (cols: SNP A1 A2 freq b se p n)
#   - reference/eQTL_besd/             GTEx v8 cis-eQTL BESD files, one set
#                                      (.besd/.esi/.epi) per tissue
#   - reference/1000G_EUR/EUR          1000 Genomes EUR PLINK bfile (LD ref)
# =============================================================================

source("config.R")
library(data.table)

GWAS_MA_DIR <- file.path(DATA_DIR, "ma")            # <TRAIT>.ma files
EQTL_DIR    <- file.path(REF_DIR, "eQTL_besd_V8")      # GTEx v8 BESD files
LD_BFILE    <- file.path(REF_DIR, "1000G_EUR", "EUR")
out_root <- file.path(RESULTS_DIR, "smr")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

all_traits <- c(MDD_TRAITS, CVD_TRAITS)

# Define target genes for SMR analysis.
# Pleiotropic SNPs supported by the colocalization analysis were functionally
# annotated, and genes located within the candidate loci were selected for
# downstream SMR analysis.
target_genes <- c("APEH","BSN","CELSR3","IP6K2","LINC03003","LOC105372666","LOC105375005","LOC124904111","MST1",
                  "PFKFB4","RHOA","RNF123","STIMATE","STIMATE-MUSTN1","TCTA","TMEM106B","TRAIP","ZNF652")
target_file <- file.path(out_root,"target_genes.txt")

write.table(target_genes,target_file,row.names = FALSE,col.names = FALSE,quote = FALSE)

# each BESD set is a prefix shared by .besd / .esi / .epi files
eqtl_prefixes <- unique(tools::file_path_sans_ext(
  list.files(EQTL_DIR, pattern = "\\.(besd|esi|epi)$", full.names = TRUE)))

## ---- Step 1: run SMR-multi + HEIDI for target genes across trait x tissue ---
for (trait in all_traits) {
  dir.create(file.path(out_root, trait), showWarnings = FALSE)
  for (eqtl in eqtl_prefixes) {
    tissue <- basename(eqtl)
    system2(SMR_BIN, c(
      "--bfile",          LD_BFILE,
      "--gwas-summary",   file.path(GWAS_MA_DIR, paste0(trait, ".ma")),
      "--beqtl-summary",  eqtl,
      "--extract-probe",  target_file,
      "--smr-multi",                     # multi-SNP-based SMR test
      "--cis-wind",        2000,          # cis window, kb
      "--peqtl-smr",      "5e-8",
      "--ld-multi-snp",   0.1,
      "--heidi-mtd",      1,
      "--diff-freq",      0.2,
      "--diff-freq-prop", 0.05,
      "--maf",            0.01,
      "--thread-num",     8,
      "--out", file.path(out_root, trait, tissue)
    ))
  }
}

## ---- Step 2: merge .msmr output across tissues, per trait -------------------
for (trait in all_traits) {
  files <- list.files(file.path(out_root, trait),
                      pattern = "\\.msmr$", full.names = TRUE)
  if (length(files) == 0) next
  merged <- data.table::rbindlist(lapply(files, function(f) {
    d <- data.table::fread(f)
    d[, tissue := tools::file_path_sans_ext(basename(f))]
    d
  }), fill = TRUE)
  data.table::fwrite(merged,
                     file.path(out_root, paste0(trait, "_merged.msmr")), sep = "\t")
}

## ---- Step 3: identify significant gene-trait associations -------------------
# Significant gene-trait associations were defined as:
# (1) FDR-adjusted SMR multi-SNP P value < 0.05
# (2) HEIDI P value > 0.05, indicating no evidence of linkage.
msmr_files <- list.files(out_root, pattern = "_merged\\.msmr$",
                         full.names = TRUE)
all_smr <- data.table::rbindlist(lapply(msmr_files, function(f) {
  d <- data.table::fread(f)
  d[, trait := sub("_merged\\.msmr$", "", basename(f))]
  d
}), fill = TRUE)


all_smr[, p_SMR_multi_FDR := p.adjust(p_SMR_multi,method = "fdr")]

sig <- all_smr[p_SMR_multi_FDR < 0.05 & p_HEIDI > 0.05]
data.table::fwrite(sig, file.path(out_root, "SMR_significant_genes.csv"))

cat("SMR + HEIDI complete.",
    nrow(sig), "significant gene-trait associations.\n")
cat("Result:", file.path(out_root, "SMR_significant_genes.csv"), "\n")
