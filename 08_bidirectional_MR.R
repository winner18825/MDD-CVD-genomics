# =============================================================================
# 08 | Bidirectional two-sample Mendelian randomization
# =============================================================================
# Bidirectional two-sample MR between the three MDD phenotypes and the
# seven cardiovascular traits.
#
# Workflow:
#   Part 1  instrument preparation -- genome-wide significant SNPs, LD
#           clumping, F-statistic filter (F > 10)
#   Part 2  per-pair MR -- harmonise -> MR-PRESSO outlier removal -> five MR
#           methods + Cochran's Q heterogeneity + MR-Egger intercept
#   Part 3  loop over all pairs in both directions, then collate results
#
# Public tools only -- no in-house wrapper required:
#   install.packages(c("TwoSampleMR", "MRPRESSO", "ieugwasr",
#                       "data.table", "dplyr"))
# Local LD clumping uses ieugwasr::ld_clump() with a PLINK binary and a
# 1000 Genomes EUR PLINK fileset (paths set in config.R).
# =============================================================================

source("config.R")
library(TwoSampleMR)
library(MRPRESSO)
library(ieugwasr)
library(data.table)

## ---- Parameters -------------------------------------------------------------
P_TH       <- 5e-6     # instrument p-value threshold
CLUMP_KB   <- 10000    # LD clumping window (kb)
CLUMP_R2   <- 0.001    # LD clumping r^2 threshold
F_FILTER   <- 10       # minimum instrument F-statistic
N_DIST     <- 1000     # MR-PRESSO simulated distributions

LD_BFILE <- file.path(REF_DIR, "1000G_EUR", "EUR")   # PLINK fileset prefix

mr_dir   <- file.path(RESULTS_DIR, "mr")
inst_dir <- file.path(mr_dir, "instruments")
dir.create(inst_dir, recursive = TRUE, showWarnings = FALSE)

all_traits <- unique(c(MDD_TRAITS, CVD_TRAITS))

## ---- Local LD clumping (public replacement for clump_data_local) -----------
clump_local <- function(dat) {
  pruned <- ieugwasr::ld_clump(
    dat        = data.frame(rsid = dat$SNP,
                            pval = dat$pval.exposure,
                            id   = dat$id.exposure),
    clump_kb   = CLUMP_KB,
    clump_r2   = CLUMP_R2,
    clump_p    = 1,
    bfile      = LD_BFILE,
    plink_bin  = PLINK_BIN)
  dat[dat$SNP %in% pruned$rsid, ]
}

# =============================================================================
# PART 1 | Instrument preparation
# =============================================================================
message("\n========== Part 1: instrument preparation ==========\n")

for (trait in all_traits) {
  message(">>> exposure: ", trait)
  input_path <- file.path(DATA_DIR, paste0(trait, ".txt"))
  if (!file.exists(input_path)) { message("  [skip] file not found");  next }
  
  exposure_dat <- tryCatch(
    read_exposure_data(
      filename = input_path, sep = "\t",
      snp_col = "SNP", beta_col = "beta", se_col = "se",
      effect_allele_col = "effect_allele", other_allele_col = "other_allele",
      eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize",
      chr_col = "chr", pos_col = "pos", clump = FALSE),
    error = function(e) { message("  [error] ", e$message); NULL })
  if (is.null(exposure_dat)) next
  # genome-wide significant SNPs
  sig <- subset(exposure_dat, pval.exposure < P_TH)
  if (nrow(sig) == 0) { message("  [skip] no significant SNP"); next }
  message("  significant SNPs: ", nrow(sig))
  
  # LD clumping (independent instruments)
  clumped <- tryCatch(clump_local(sig),
                      error = function(e) { message("  [error] clumping: ", e$message); NULL })
  if (is.null(clumped) || nrow(clumped) == 0) next
  message("  after clumping: ", nrow(clumped))
  
  # instrument strength: per-SNP R^2 and F-statistic, keep F > F_FILTER
  clumped <- transform(clumped,
                       R2 = 2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure))
  clumped <- transform(clumped,
                       F  = (samplesize.exposure - 2) * R2 / (1 - R2))
  outTab <- clumped[clumped$F > F_FILTER, ]
  if (nrow(outTab) == 0) { message("  [skip] no SNP after F filter"); next }
  message("  after F filter: ", nrow(outTab))
  
  fwrite(outTab,
         file.path(inst_dir, paste0(trait, ".instruments.txt")), sep = "\t")
}

# =============================================================================
# PART 2 | MR analysis for one exposure -> outcome pair
# =============================================================================
run_one_MR <- function(trait1, trait2) {
  message("\n------ MR: ", trait1, " -> ", trait2, " ------")
  
  exposureFile <- file.path(inst_dir, paste0(trait1, ".instruments.txt"))
  outcome_path <- file.path(DATA_DIR, paste0(trait2, ".txt"))
  if (!file.exists(exposureFile) || !file.exists(outcome_path)) {
    message("  [skip] input file not found"); return(NULL)
  }
  
  pair_dir <- file.path(mr_dir, paste0(trait1, "_", trait2))
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- read exposure / outcome ----------------------------------------------
  exposure_dat <- read_exposure_data(
    filename = exposureFile, sep = "\t",
    snp_col = "SNP",
    beta_col = "beta.exposure", se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure", pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure",
    chr_col = "chr.exposure", pos_col = "pos.exposure", clump = FALSE)
  if (nrow(exposure_dat) == 0) return(NULL)
  
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP, filename = outcome_path, sep = "\t",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize",
    chr_col = "chr", pos_col = "pos")
  if (nrow(outcome_dat) == 0) return(NULL)
  
  exposure_dat$exposure <- exposure_dat$id.exposure <- trait1
  outcome_dat$outcome   <- outcome_dat$id.outcome   <- trait2
  
  # ---- harmonise ------------------------------------------------------------
  dat_raw <- harmonise_data(exposure_dat, outcome_dat)
  if (nrow(dat_raw) == 0) return(NULL)
  message("  SNPs after harmonisation: ", nrow(dat_raw))
  
  # ---- Step 1: MR-PRESSO outlier detection / removal ------------------------
  presso_result <- NULL
  outlier_snps  <- character(0)
  req <- c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
  dat <- dat_raw[complete.cases(dat_raw[, req]), ]
  
  if (nrow(dat) >= 4) {
    presso_result <- tryCatch({
      set.seed(1234)
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",     SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = as.data.frame(dat),
        NbDistribution = N_DIST, SignifThreshold = 0.05)
    }, error = function(e) { message("  [warn] MR-PRESSO: ", e$message); NULL })
    
    if (!is.null(presso_result)) {
      save(presso_result, file = file.path(pair_dir,
                                           paste0(trait1, "_", trait2, ".MR-PRESSO_result.RData")))
      ot <- presso_result$`MR-PRESSO results`$`Outlier Test`
      if (!is.null(ot) && nrow(ot) > 0) {
        fwrite(ot, file.path(pair_dir,
                             paste0(trait1, "_", trait2, ".MR-PRESSO_outliers.txt")),
               row.names = TRUE, sep = "\t")
        # flag outliers at a Bonferroni-corrected threshold
        idx <- which(ot$Pvalue < 0.05 / nrow(dat))
        if (length(idx) > 0) {
          outlier_snps <- dat$SNP[idx]
          dat <- dat[!dat$SNP %in% outlier_snps, ]
          message("  MR-PRESSO removed ", length(outlier_snps), " outlier SNP(s)")
        }
      }
    }
  } else {
    message("  [note] < 4 SNPs, MR-PRESSO skipped")
  }
  message("  SNPs after outlier removal: ", nrow(dat))
  fwrite(dat, file.path(pair_dir,
                        paste0(trait1, "_", trait2, ".instruments_final.txt")), sep = "\t")
  if (nrow(dat) < 3) return(list(mr_res = NULL, het_res = NULL,
                                 plt_res = NULL, presso = presso_result,
                                 outlier_snps = outlier_snps))
  
  # ---- Step 2: five MR methods ----------------------------------------------
  mr_res <- mr(dat, method_list = c(
    "mr_ivw", "mr_egger_regression", "mr_weighted_median",
    "mr_weighted_mode", "mr_simple_mode"))
  if (!is.null(mr_res) && nrow(mr_res) > 0) {
    mr_res$direction <- paste0(trait1, " -> ", trait2)
    mr_res$OR        <- exp(mr_res$b)
    mr_res$OR_lci95  <- exp(mr_res$b - 1.96 * mr_res$se)
    mr_res$OR_uci95  <- exp(mr_res$b + 1.96 * mr_res$se)
    fwrite(mr_res, file.path(pair_dir,
                             paste0(trait1, "_", trait2, ".MR_5methods.txt")), sep = "\t")
  }
  
  # ---- Step 3: heterogeneity (Cochran's Q) ----------------------------------
  het_res <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  if (!is.null(het_res) && nrow(het_res) > 0) {
    het_res$direction <- paste0(trait1, " -> ", trait2)
    fwrite(het_res, file.path(pair_dir,
                              paste0(trait1, "_", trait2, ".heterogeneity.txt")), sep = "\t")
  }
  
  # ---- Step 4: directional pleiotropy (MR-Egger intercept) ------------------
  plt_res <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if (!is.null(plt_res) && nrow(plt_res) > 0) {
    plt_res$direction <- paste0(trait1, " -> ", trait2)
    fwrite(plt_res, file.path(pair_dir,
                              paste0(trait1, "_", trait2, ".pleiotropy.txt")), sep = "\t")
  }
  
  message("  [done] ", trait1, " -> ", trait2)
  list(mr_res = mr_res, het_res = het_res, plt_res = plt_res,
       presso = presso_result, outlier_snps = outlier_snps)
}

# =============================================================================
# PART 3 | Run every pair in both directions and collate
# =============================================================================
message("\n========== Part 3: bidirectional MR ==========\n")

all_mr <- list(); all_het <- list(); all_plt <- list(); all_presso <- list()
log_rows <- list()

collect <- function(t1, t2) {
  res <- tryCatch(run_one_MR(t1, t2),
                  error = function(e) { message("  [fatal] ", e$message); NULL })
  log_rows[[length(log_rows) + 1]] <<- data.frame(
    exposure = t1, outcome = t2,
    status = if (is.null(res)) "FAILED" else "SUCCESS")
  if (is.null(res)) return(invisible())
  if (!is.null(res$mr_res))  all_mr[[length(all_mr) + 1]]   <<- res$mr_res
  if (!is.null(res$het_res)) all_het[[length(all_het) + 1]] <<- res$het_res
  if (!is.null(res$plt_res)) all_plt[[length(all_plt) + 1]] <<- res$plt_res
  gt <- res$presso$`MR-PRESSO results`$`Global Test`
  if (!is.null(gt)) {
    pv <- suppressWarnings(as.numeric(gt$Pvalue))
    if (is.na(pv)) pv <- 1 / N_DIST          # below-resolution global-test p
    all_presso[[length(all_presso) + 1]] <<- data.frame(
      direction = paste0(t1, " -> ", t2),
      PRESSO_RSSobs = gt$RSSobs, PRESSO_Global_pval = pv,
      N_outliers = length(res$outlier_snps))
  }
}

for (md in MDD_TRAITS) {
  for (cv in CVD_TRAITS) {
    collect(md, cv)   # MDD -> CVD
    collect(cv, md)   # CVD -> MDD
  }
}

## ---- Collated tables --------------------------------------------------------
fwrite(rbindlist(log_rows, fill = TRUE),
       file.path(mr_dir, "MR_analysis_log.txt"), sep = "\t")

if (length(all_mr) > 0) {
  df_all <- rbindlist(all_mr, fill = TRUE)
  fwrite(df_all, file.path(mr_dir, "MR_all_methods_results.txt"), sep = "\t")
  
  # IVW summary 
  df_ivw <- df_all[method == "Inverse variance weighted"]

  fwrite(df_ivw, file.path(mr_dir, "MR_IVW_results_summary.txt"), sep = "\t")
}
if (length(all_het) > 0)
  fwrite(rbindlist(all_het, fill = TRUE),
         file.path(mr_dir, "MR_heterogeneity_results.txt"), sep = "\t")
if (length(all_plt) > 0)
  fwrite(rbindlist(all_plt, fill = TRUE),
         file.path(mr_dir, "MR_pleiotropy_results.txt"), sep = "\t")
if (length(all_presso) > 0)
  fwrite(rbindlist(all_presso, fill = TRUE),
         file.path(mr_dir, "MR_PRESSO_global_test_summary.txt"), sep = "\t")

cat("\nBidirectional MR complete. Results in:", mr_dir, "\n")

