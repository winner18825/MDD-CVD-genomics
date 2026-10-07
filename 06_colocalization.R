# =============================================================================
# 06 | Bayesian colocalization (coloc)
# =============================================================================
# For every pleiotropic lead SNP identified in script 05, tests whether the
# two traits share a single causal variant within a +/- 500 kb window using
# coloc.abf (Giambartolomei et al., PLoS Genet 2014) with default priors
# (p1 = p2 = 1e-4, p12 = 1e-5). Colocalization is declared when PP.H4 > 0.7.
#
# coloc.abf forms a per-SNP approximate Bayes factor; SNPs are matched across
# the two GWAS by rsID. Public package only:  install.packages("coloc")
#
# This script also produces the consolidated summary table (replacing the
# separate per-locus file-merging step of the original PPH.R).
# =============================================================================

source("config.R")
library(data.table)
library(coloc)

WINDOW        <- 500000   # +/- bp around each lead SNP
PP4_THRESHOLD <- 0.8      # posterior probability cut-off for colocalization

# Case proportion (n_cases / n_total) for each trait -- required for type "cc".
# Replace with the values for the GWAS datasets you use.
CASE_PROP <- c(
  MDD_ALL              = 0.294646,
  eoMDD                = 0.304223,
  loMDD                = 0.234368,
  Hypertension         = 0.288834,
  AF                   = 0.058807,
  AIS                  = 0.077708,
  CAD                  = 0.155721,
  HF                   = 0.081677,
  MI                   = 0.045292,
  PAD                  = 0.017905)

out_dir <- file.path(RESULTS_DIR, "coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

loci <- data.table::fread(
  file.path(RESULTS_DIR, "merge_cpassoc_placo", "merged_all_placo_sig.txt"))

# read each full GWAS only once
gwas_cache <- new.env()
get_gwas <- function(trait) {
  if (!exists(trait, envir = gwas_cache))
    assign(trait,
           data.table::fread(file.path(DATA_DIR, paste0(trait, ".standard.txt"))),
           envir = gwas_cache)
  get(trait, envir = gwas_cache)
}

make_dataset <- function(df, trait) {
  list(snp     = df[[COL$SNP]],
       beta    = df[[COL$BETA]],
       varbeta = df[[COL$SE]]^2,
       MAF     = df[[COL$EAF]],
       N       = df[[COL$N]],
       type    = "cc",
       s       = CASE_PROP[[trait]])
}

summary_rows <- list()
for (i in seq_len(nrow(loci))) {
  pair <- strsplit(loci$prefix[i], "&", fixed = TRUE)[[1]]
  t1 <- pair[1]; t2 <- pair[2]
  snp <- loci$SNP[i]
  chr <- loci[[paste0(COL$CHR, ".x")]][i]
  bp  <- loci[[paste0(COL$BP,  ".x")]][i]
  
  g1 <- get_gwas(t1); g2 <- get_gwas(t2)
  reg1 <- g1[g1[[COL$CHR]] == chr & abs(g1[[COL$BP]] - bp) <= WINDOW, ]
  reg2 <- g2[g2[[COL$CHR]] == chr & abs(g2[[COL$BP]] - bp) <= WINDOW, ]
  
  common <- intersect(reg1[[COL$SNP]], reg2[[COL$SNP]])
  if (length(common) < 50) next                 # too few SNPs to colocalize
  reg1 <- reg1[match(common, reg1[[COL$SNP]]), ]
  reg2 <- reg2[match(common, reg2[[COL$SNP]]), ]

  res <- coloc::coloc.abf(
    dataset1 = make_dataset(reg1, t1),
    dataset2 = make_dataset(reg2, t2),
    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  
  s <- res$summary
  summary_rows[[length(summary_rows) + 1]] <- data.table::data.table(
    trait1 = t1, trait2 = t2, lead_SNP = snp, CHR = chr, BP = bp,
    nsnps = s[["nsnps"]],
    PP.H0 = s[["PP.H0.abf"]], PP.H1 = s[["PP.H1.abf"]],
    PP.H2 = s[["PP.H2.abf"]], PP.H3 = s[["PP.H3.abf"]],
    PP.H4 = s[["PP.H4.abf"]])
}

coloc_summary <- data.table::rbindlist(summary_rows)
coloc_summary[, colocalized := PP.H4 > PP4_THRESHOLD]
data.table::fwrite(coloc_summary,
                   file.path(out_dir, "coloc_summary.txt"), sep = "\t")

cat("Colocalization complete.",
    sum(coloc_summary$colocalized), "of", nrow(coloc_summary),
    "loci with PP.H4 >", PP4_THRESHOLD, "\n")
cat("Summary table:", file.path(out_dir, "coloc_summary.txt"), "\n")
