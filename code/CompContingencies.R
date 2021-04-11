CompContingencies <- function(tumor_mat, normal_mat,
                              includeORs = F){
  # Uses a fisher test based framework to compute 2 2x2 fisher tests
  # and then compares these 2 fisher tests results
  # Assumes normality, which is a simplifying assumption
  
  # if include ORs is T, include the tumor matrix OR,
  # and the normal matrix OR in the output vector
  # as 3rd and 4th items
  
  # read packages and or other functions
  library(BSDA)
  source("SElogodds.R")
  
  # run fishers exact test to get odds ratios and confidence intervals for each
  # contigency table matrix
  tumor_fish <- fisher.test(tumor_mat)
  normal_fish <- fisher.test(normal_mat)
  
  tumor_lodds <- log(tumor_fish$estimate)
  normal_lodds <- log(normal_fish$estimate)
  tumor_SE_logodds <- SElogodds(tumor_mat)
  normal_SE_logodds <- SElogodds(normal_mat)
  delta <- tumor_lodds - normal_lodds
  names(delta) <- NULL
  SE_delta <- sqrt(tumor_SE_logodds**2 + normal_SE_logodds**2)
  
  z <- delta / SE_delta
  pval <- pnorm(-abs(z))
  
  if (isTRUE(includeORs)){
    return(c(delta, pval,tumor_fish$estimate,
             normal_fish$estimate))
  } else{
    return(c(delta, pval))
  }
}
