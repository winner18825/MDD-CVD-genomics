# =============================================================================
# 03 | Cross-trait meta-analysis -- CPASSOC (SHet)
# =============================================================================
# Detects pleiotropic loci shared by a MDD phenotype and a CVD trait
# using the SHet statistic from CPASSOC (Li et al., Methods Mol Biol 2017).
# The SHet statistic accounts for cross-trait association signals by jointly
# analyzing summary statistics and is designed to accommodate heterogeneity
# of effect sizes across phenotypes.
#
# The CPASSOC functions used in this analysis were implemented based on the
# original method description and workflow. The required functions are provided
# locally in: CPASSOC_FunctionSet.R
#
# This implementation enables reproducible application of the CPASSOC SHet
# statistic without requiring external package installation.
# =============================================================================

source("config.R")
source("CPASSOC_FunctionSet.R")
library(dplyr)
library(stringr)
library(data.table)
library(future)

out_dir <- file.path(RESULTS_DIR, "CPASSOC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_cpassoc <- function(mdd, cvd){
  trait1_path <- file.path(DATA_DIR,paste0(mdd,".txt"))
  trait2_path <- file.path(DATA_DIR,paste0(cvd,".txt"))
  
  trait1_data <- fread(trait1_path, data.table = FALSE)
  trait2_data <- fread(trait2_path, data.table = FALSE)
  
  trait1_samplesize <- max(trait1_data$N, na.rm = TRUE)
  trait2_samplesize <- max(trait2_data$N, na.rm = TRUE)
  
  sample_size <- c(trait1_samplesize, trait2_samplesize)
  m <- cpassoc.combine_assoc(files_path = c(trait1_path,trait2_path),traits = c(mdd,cvd))
  # Filter SNPs: select target SNPs and filter based on Z-scores
  cor_mat <- m %>%
    # Filter 1: retain independent SNPs
    dplyr::filter(SNP %in%  snplist$V1) %>%
    # Filter 2: remove SNPs with |Z-score| > 1.96
    dplyr::filter(dplyr::if_all(tidyselect::starts_with('z'), ~ abs(.x) < 1.96)) %>%
    dplyr::select(tidyselect::starts_with('z')) %>%
    as.matrix() %>%
    cor()
  # CPASSOC
  sh <- cpassoc.shet_shom_test(gwas_data = m,SHet_test = TRUE,SHom_test = TRUE,workers = 4,mvrnorm_n = 1e6,set_seed = 666)
  # Save results
  fwrite(sh,file = paste0(out_dir,"/",mdd,"_",cvd,"_CPASSOC_SHet.csv"),quote = FALSE,sep="\t")
}

# plink 
sig_results_all <- list()
# Use PLINK to select a set of independent SNPs
# The generated list of independent SNPs can be reused for the same analysis scenario
args <- c("--bfile",file.path(PLINK_BIN,"1kg.v3/EUR"),
          "--indep-pairwise",500,50,0.2,
          "--out","plink_pruning")
system2(command = plink,args = args)

# The above code generates two output files
# plink_pruning.prune.in   # Independent SNPs retained after pruning
# plink_pruning.prune.out  # SNPs removed during pruning
snplist <-  data.table::fread(input = paste0("plink_pruning",".prune.in"),header = F,data.table = F)

for(mdd in MDD_TRAITS){
  for(cvd in CVD_TRAITS){
    result <- run_cpassoc(mdd,cvd)
  }
}

cat("CPASSOC cross-trait meta-analysis complete. Results in:", out_dir, "\n")