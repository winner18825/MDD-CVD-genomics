# =============================================================================
# 05 | Merge CPASSOC + PLACO significant loci
# =============================================================================
# Genome-wide significant loci identified by CPASSOC and PLACO were integrated
# to obtain variants consistently supported by both complementary pleiotropy
# detection approaches.

# The overlap between CPASSOC and PLACO significant SNP sets was obtained,
# followed by LD clumping using PLINK with 1000 Genomes European reference
# variants (r2 < 0.2, window size 500 kb). The resulting independent
# pleiotropic loci were combined across all trait pairs and used for
# downstream colocalization analyses.
# =============================================================================

source("config.R")
library(data.table)

cpassoc_dir   <- file.path(RESULTS_DIR, "CPASSOC")
placo_dir     <- file.path(RESULTS_DIR, "placo")
out_dir       <- file.path(RESULTS_DIR, "merge_cpassoc_placo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged_all <- list()
for (m in MDD_TRAITS) {
  for (cv in CVD_TRAITS) {
    tryCatch({
    cp_file <- file.path(cpassoc_dir, paste0(m, "_", cv, "_cpassoc_result.txt"))
    pl_file <- file.path(placo_dir,   paste0(m, "_", cv, "_PLACO_result.txt"))
    if (!file.exists(cp_file) || !file.exists(pl_file)) next
      
    cp <- data.table::fread(cp_file)
    pl <- data.table::fread(pl_file)
      
    # CPASSOC: genome-wide significant SHet, both single-trait p < 1e-3
    cp <- cp[p_Shet < 5e-8 & pval.x < 1e-3 & pval.y < 1e-3]
    # PLACO: genome-wide significant composite p
    pl <- pl[p.placo < 5e-8]
    
    shared <- merge(cp, pl, by = "SNP", suffixes = c(".x", ".y"))
    # plink 
    fn <- tempfile()
    write.table(data.frame(SNP=shared[["SNP"]], P=shared[["p_SHet"]]), file=fn, row.names=F, col.names=T, quote=F)
    args <- c("--bfile",file.path(PLINK_BIN,"1kg.v3/EUR"),
              "--clump",fn,
              "--clump-p1",5e-8,
              "--clump-p2",1e-5,
              "--clump-r2",0.2,
              "--clump-kb",500,
              "--out","plink_clump")
    system2(command = plink,args = args)
    unlink(fn)
    
    snplist_clump <- data.table::fread("plink_clump.clumped",data.table = F)
    snplist_clump <- snplist_clump$SNP
    
    shared <- shared[shared$SNP %in% snplist_clump, ]
    if (nrow(shared) == 0) next
    
    shared[, prefix := paste0(m, "&", cv)]
    data.table::fwrite(shared,
                       file.path(out_dir, paste0(m, "_", cv, "_merge_cpassoc_placo_sig.txt")),
                       sep = "\t")
    merged_all[[paste0(m, "_", cv)]] <- shared
    }, error=function(e){
      message("ERROR: ",m," - ",cv," : ",e$message)
    })
  }
}

all_loci <- data.table::rbindlist(merged_all, fill = TRUE)
data.table::fwrite(all_loci,
                   file.path(out_dir, "merged_all_placo_sig.txt"), sep = "\t")

cat("Merged", nrow(all_loci), "pleiotropic loci across",
    length(merged_all), "trait pairs.\n")
cat("Combined table:",
    file.path(out_dir, "merged_all_placo_sig.txt"), "\n")


  