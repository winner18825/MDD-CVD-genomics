# =============================================================================
# 02 | Global genetic correlation -- HDL
# =============================================================================
# High-Definition Likelihood (HDL; Ning, Pawitan & Shen, Nat Genet 2020)
# estimate of genome-wide genetic correlation, used as a complementary,
# higher-precision method alongside LDSC.
#
# HDL is a public R package -- no in-house wrapper required:
#   remotes::install_github("zhenin/HDL/HDL")
# The EUR LD reference panel (UK Biobank imputed, SVD eigen99 extraction) is
# distributed with HDL:  https://github.com/zhenin/HDL
# =============================================================================

source("config.R")
library(HDL)
library(data.table)

# HDL LD reference panel
HDL_LD_PATH <- file.path(REF_DIR, "HDL", "UKB_imputed_SVD_eigen99_extraction")

out_dir <- file.path(RESULTS_DIR, "hdl")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# HDL requires columns SNP, A1, A2, N, Z. Loads <TRAIT>_MTAG.txt and
# harmonises common header variants.
load_for_hdl <- function(trait) {
  df <- as.data.frame(data.table::fread(
    file.path(DATA_DIR, paste0(trait, "_MTAG.txt"))))
  ren <- c(effect_allele = "A1", other_allele = "A2")
  hit <- names(df) %in% names(ren)
  names(df)[hit] <- ren[names(df)[hit]]
  if (!"Z" %in% names(df) && all(c("beta", "se") %in% names(df)))
    df$Z <- df$beta / df$se
  df[, c("SNP", "A1", "A2", "N", "Z")]
}

results <- list()
for (m in MDD_TRAITS) {
  g1 <- load_for_hdl(m)
  for (cv in CVD_TRAITS) {
    g2  <- load_for_hdl(cv)
    res <- HDL.rg(
      gwas1.df    = g1,
      gwas2.df    = g2,
      LD.path     = HDL_LD_PATH,
      Nref        = 335265,
      eigen.cut   = "automatic",
      output.file = file.path(out_dir, paste0(m, "_", cv, "_HDL.log"))
    )
    results[[paste0(m, "_", cv)]] <- data.frame(
      MDD = m, cvd = cv,
      rg = res$rg, rg_se = res$rg.se, P = res$P)
  }
}

rg_table <- do.call(rbind, results)
data.table::fwrite(rg_table, file.path(out_dir, "HDL_rg_summary.csv"))
cat("HDL genetic correlation complete. Summary:",
    file.path(out_dir, "HDL_rg_summary.csv"), "\n")