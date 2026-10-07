# =============================================================================

library("MASS")
library("Matrix")
library("compiler")

Non_Trucated_TestScore <- function(X, SampleSize, CorrMatrix) {
  Wi = matrix(SampleSize, nrow = 1);
  sumW = sqrt(sum(Wi^2));
  W = Wi / sumW;

  Sigma = ginv(CorrMatrix);
  XX = apply(X, 1, function(x) {
    x1 <- matrix(x, ncol = length(x), nrow = 1);
    Tmat = W %*% Sigma %*% t(x1);
    Tmat = (Tmat*Tmat) / (W %*% Sigma %*% t(W));
    return(Tmat[1,1]);
  }
  );
  return(XX);
}
SHom <- cmpfun(Non_Trucated_TestScore);

Trucated_TestScore <- function(
    X, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05,
    isAllpossible = T
) {
  N <- nrow(X)

  Wi <- matrix(SampleSize, nrow = 1)
  sumW <- sqrt(sum(Wi^2))
  W <- Wi / sumW

  calc_TTT <- function(x) {
    TTT = -1

    if (isAllpossible == T ) {
      cutoff <- sort(unique(abs(x)))	  ## it will filter out any of them.
    } else {
      cutoff = seq(startCutoff, endCutoff, CutoffStep)
    }

    for (threshold in cutoff) {
      x1 <- x
      index <- which(abs(x1) < threshold)

      if (length(index) == N) break;

      A <- CorrMatrix

      W1 <- W
      if (length(index) != 0) {
        x1 <- x1[-index]
        A  <- A[-index, -index]   ## update the matrix
        W1 <- W[-index]
      }

      if (correct == 1)	{
        W1 <- W1 * sign(x1)
      }

      A <- ginv(A)
      x1 <- matrix(x1, nrow = 1)
      W1 <- matrix(W1, nrow = 1)
      Tmat <- W1 %*% A %*% t(x1)
      Tmat <- (Tmat * Tmat) / (W1 %*% A %*% t(W1))
 
      if (is.na(Tmat[1, 1])) {
        TTT <- NA
      } else {
        if (TTT < Tmat[1, 1]) {
          TTT <- Tmat[1, 1]
        }
      }
    }
    return(TTT)
  }
  res <- future.apply::future_apply(X, 1, calc_TTT, future.seed = TRUE)

  return(res)
}
SHet <- cmpfun(Trucated_TestScore);

EstimateGamma <- function (N = 1E6, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = T) {

  Wi = matrix(SampleSize, nrow = 1);
  sumW = sqrt(sum(Wi^2));
  W = Wi / sumW;

  Permutation = mvrnorm(n = N, mu = c(rep(0, length(SampleSize))), Sigma = CorrMatrix, tol = 1e-8, empirical = F);

  Stat =  Trucated_TestScore(X = Permutation, SampleSize = SampleSize, CorrMatrix = CorrMatrix,
                             correct = correct, startCutoff = startCutoff, endCutoff = endCutoff,
                             CutoffStep = CutoffStep, isAllpossible = isAllpossible);
  a = min(Stat)*3/4
  ex3 = mean(Stat*Stat*Stat)
  V =	var(Stat);

  for (i in 1:100){
    E = mean(Stat)-a;
    k = E^2/V
    theta = V/E
    a = (-3*k*(k+1)*theta**2+sqrt(9*k**2*(k+1)**2*theta**4-12*k*theta*(k*(k+1)*(k+2)*theta**3-ex3)))/6/k/theta
  }

  para = c(k,theta,a);
  return(para);
}

######## Combine GWAS summary statistics
cpassoc.combine_assoc <- function(files_path,traits=NULL) {
  
  if(!is.null(traits) & (length(traits) == length(files_path)))
  {
    traits <- traits
  }else{
    traits <- NULL
  }
  
  gwas_data <- data.frame()
  
  for (index in 1:length(files_path)) {
    dat <- data.table::fread(files_path[index], data.table = F, fill = T)
    
    # Filter and clean SNP identifiers (rs IDs)
    dat <- dat %>% dplyr::filter(grepl("^rs", SNP))
    dat <- dat %>% dplyr::mutate(SNP = sub(",.*$", "", SNP))
    dat <- dat %>% dplyr::distinct(SNP, .keep_all = TRUE)
    # Calculate Z-scores
    dat <-  dat %>% dplyr::mutate(z = beta / se)
    
    col <- c("SNP", "chr", "pos", "effect_allele", "other_allele","z")
    dat <- dat %>%  dplyr::select(all_of(col))
    
    new_names <-  c("SNP","chr","pos","effect_allele", "other_allele")
    
    if(!is.null(traits)){
      z <- paste0("z_",traits[index])
    }else{
      z <- paste0("z_", index)
    }
    
    colnames(dat) <- c(new_names, z)
    if (index == 1) {
      gwas_data <- dat
    } else{
      
      col <- c("SNP", "effect_allele", "other_allele", z)
      dat <- dat %>%  dplyr::select(col)
      dat <- dat %>%  dplyr::filter(SNP %in% gwas_data$SNP)
      effect_allele <- paste0("effect_allele", index)
      other_allele <- paste0("other_allele", index)
      colnames(dat) <- c("SNP", effect_allele, other_allele,z )
      
      gwas_data <- dplyr::inner_join(gwas_data, dat, by = "SNP")
      
      # Resolve inconsistent effect_allele and other_allele orientations
      snps_to_flip <- (gwas_data$effect_allele == gwas_data[, other_allele]) &
        (gwas_data$other_allele == gwas_data[, effect_allele])
      snps_to_keep <- (gwas_data$effect_allele == gwas_data[, effect_allele]) &
        (gwas_data$other_allele == gwas_data[, other_allele])
      snps_to_keep <- snps_to_keep | snps_to_flip
      
      gwas_data <- gwas_data[snps_to_keep, ]
      # After merging, check again for reversed allele orientations
      snps_to_flip <- (gwas_data$effect_allele == gwas_data[, other_allele]) &
        (gwas_data$other_allele == gwas_data[, effect_allele])
      if (sum(snps_to_flip) > 0) {
        # Swap the data positions
        diff_index <- which(snps_to_flip)
        gwas_data[diff_index, z] <- -1 * gwas_data[diff_index , z]
      }
      
      # Remove unnecessary columns
      gwas_data[, effect_allele] <- NULL
      gwas_data[, other_allele] <- NULL
    }
  }
  
  return (gwas_data)
}

########## shet_shom_test
cpassoc.shet_shom_test <- function(gwas_data,SHet_test =T,SHom_test =T,workers =4,set_seed = 666,mvrnorm_n =1e6){
  # Select Z-score columns
  Zscores <- dplyr::select(d, starts_with('z'))
  
  # SHet test
  if (SHet_test) {
    cat('Estimate parameters of gamma distribution...\n')
    
    set.seed(set_seed)
    para <- EstimateGamma(N = mvrnorm_n,
                          SampleSize = sample_size,
                          CorrMatrix = cor_mat)
    
    cat('Calculating SHet statistics...\n')
    old_plan <- future::plan(strategy = "multisession", workers = workers)
    x <- SHet(X = Zscores, SampleSize = sample_size, CorrMatrix = cor_mat, correct = 1, isAllpossible = T)
    future::plan(old_plan)
    
    cat('Calculating p-values...\n')
    d$p_SHet <- pgamma(q = x - para[3], shape = para[1], scale = para[2], lower.tail = F)
  }
  
  # SHom_test test
  if (SHom_test) {
    
    cat('Calculating SHom statistics...\n')
    SHom_res <- SHom(X = Zscores, SampleSize = sample_size, CorrMatrix = cor_mat)
    
    d$p_SHom <- pchisq(SHom_res, df = 1, ncp = 0, lower.tail = F)
    
  }
  
  # Return results
  return(d)
}





