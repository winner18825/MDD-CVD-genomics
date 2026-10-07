# =============================================================================
# 04 | Cross-trait meta-analysis -- PLACO
# =============================================================================
# Detects pleiotropic loci under a composite null hypothesis using PLACO
# (Ray & Chatterjee, PLoS Genet 2020).
#
# PLACO is distributed as a single R source file. Download the official
# version and place it at  tools/PLACO/PLACO.R :
#   https://github.com/RayDebashree/PLACO
#
# IMPORTANT: follow the official README. Verify the argument names of
# var.placo() and placo(), and the names of the returned statistics, against
# the version you download -- the workflow below reflects the documented
# interface but the authors may revise it.
# =============================================================================

source("config.R")
library(data.table)
library(parallel)

source(file.path("tools", "PLACO", "PLACO.R"))   # provides var.placo() / placo()

out_dir <- file.path(RESULTS_DIR, "placo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

harmonize <- function(m){
  # Determine whether the effect alleles are aligned in opposite directions
  m[, flip := 
      effect_allele.x == other_allele.y &
      other_allele.x == effect_allele.y
  ]
  m[, same :=
      effect_allele.x == effect_allele.y &
      other_allele.x == other_allele.y
  ]
  # Keep SNPs with compatible allele configurations
  m <- m[same | flip]
  # Note: the flip indicator is updated together with the filtered data table
  m[flip == TRUE, beta.y := -beta.y]
  # Remove auxiliary columns used for allele matching
  m[, c("same","flip") := NULL]
  return(m)
}

run_placo <- function(trait1, trait2) {
  d1 <- data.table::fread(file.path(DATA_DIR, paste0(trait1, "_MTAG.txt")))
  d2 <- data.table::fread(file.path(DATA_DIR, paste0(trait2, "_MTAG.txt")))
  m  <- merge(d1, d2, by = COL$SNP, suffixes = c(".x", ".y"))
  m <- harmonize(m)
  # Z-score matrix (compute from beta/se if a Z column is absent)
  if (all(c("Z.x", "Z.y") %in% names(m))) {
    Z <- cbind(m$Z.x, m$Z.y)
  } else {
    Z <- cbind(m[[paste0(COL$BETA, ".x")]] / m[[paste0(COL$SE, ".x")]],
               m[[paste0(COL$BETA, ".y")]] / m[[paste0(COL$SE, ".y")]])
  }
  # corresponding p-value matrix (derive from Z if absent)
  P <- 2 * pnorm(-abs(Z))
  
  # remove variants with extreme squared Z (PLACO recommendation)
  keep <- rowSums(Z^2 > 80) == 0
  m <- m[keep, ]; Z <- Z[keep, , drop = FALSE]; P <- P[keep, , drop = FALSE]
  
  # estimate variance parameters, then test every variant
  VarZ <- var.placo(Z, P, p.threshold = 1e-4)
  res  <- t(vapply(seq_len(nrow(Z)),
                   function(i) unlist(placo(Z = Z[i, ], VarZ = VarZ)),
                   numeric(2)))
  
  data.table::data.table(
    SNP     = m[[COL$SNP]],
    CHR     = m[[paste0(COL$CHR, ".x")]],
    BP      = m[[paste0(COL$BP, ".x")]],
    T.placo = res[, "T.placo"],
    p.placo = res[, "p.placo"])
}

pairs <- expand.grid(trait1 = MDD_TRAITS,trait2 = CVD_TRAITS,stringsAsFactors = FALSE)
mclapply(seq_len(nrow(pairs)),
  function(i){
    trait1 <- pairs$trait1[i]
    trait2 <- pairs$trait2[i]
    message("Running: ",trait1," vs ",trait2)
    res <- run_placo(trait1, trait2)
    outfile <- file.path(out_dir,paste0(trait1,"_",trait2,"_PLACO_result.txt"))
    data.table::fwrite(res,outfile,sep="\t")
    return(outfile)
  },
  mc.cores = 8
)

cat("PLACO cross-trait meta-analysis complete. Results in:", out_dir, "\n")