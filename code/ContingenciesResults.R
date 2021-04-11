ContingenciesResults <- function(vecs,includeORs = F){
  # takes a list of vectors (vecs) with 8 counts in each vector, and assumes the first four
  # counts per vector are for the first 2x2 table and the second four counts per
  # vector are for the second 2x2 table Runs my CompContingencies function and
  # stores deltas and pvals in a df. names(vec) is assumed to be TFs, which is how I wrote the counts parsing function
  
  # if include ORs is T, includes the tumor and normal
  # matrix ORs in the output df
  
  source("CompContingencies.R")
  
  # empty vectors
  TFs <- c()
  deltas <- c()
  pvals <- c()
  tumor_ORs <- c()
  normal_ORs <- c()
  
  # iterate
  for (i in 1:length(vecs)){
    mat1 <- matrix(vecs[[i]][1:4], nrow = 2)
    mat2 <- matrix(vecs[[i]][5:8], nrow = 2)
    
    res <- CompContingencies(mat1, mat2, includeORs = T)
    delta <- res[1]
    pval <- res[2]
    tumor_OR <- res[3]
    normal_OR <- res[4]
    
    TFs <- append(TFs, names(vecs)[i])
    deltas <- append(deltas, delta)
    pvals <- append(pvals, pval)
    tumor_ORs <- append(tumor_ORs, tumor_OR)
    normal_ORs <- append(normal_ORs, normal_OR)
  }
  
  # make output df
  output <- data.frame(TF = TFs,
                       delta = deltas,
                       pval = pvals,
                       tumor_OR = tumor_ORs,
                       normal_OR = normal_ORs)
  # return output df
  return(output)
}