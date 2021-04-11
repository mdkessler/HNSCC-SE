TFTargetJaccard <- function(TF_gene_list, dist = T){
  # TF_gene_list is a list, with TF as keys and DF as vals
  # The val df contains one col, with target genes of TF
  # This col is named TF_genes
  
  # This function will calculate the jaccard based overlap
  # of the target genes between all pairwise TFs
  
  # The function will return a matrix with Tfs as cols and rows
  # (i.e. a pairwise matrix) with these jaccard values
  
  # if the dist = T argument, which is default, than 1-jaccard
  # will be returned as the matrix values instead
  
  # import packages
  library(OmicsMarkeR) # has jaccard function
  output_mat <- c() # initializing
  for (TF1 in names(TF_gene_list)){
    TF1_col <- c() # initializing
    for (TF2 in names(TF_gene_list)){
      x <- as.character(TF_gene_list[[TF1]]$TF_genes)
      y <- as.character(TF_gene_list[[TF2]]$TF_genes)
      jac <- jaccard(x, y)
      if (isTRUE(dist)){
        jac <- 1 - jac
      }
      TF1_col <- append(TF1_col, jac)
    }
    # now store the TF1_col of jac/dist values in output_mat
    output_mat <- cbind(output_mat, TF1_col)
  }
  # set row and col names
  colnames(output_mat) <- names(TF_gene_list)
  rownames(output_mat) <- names(TF_gene_list)
  return(output_mat)
}