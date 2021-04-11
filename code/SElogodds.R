SElogodds <- function(mat){
  SE <- sqrt((1/mat[1,1]) + (1/mat[1,2]) + (1/mat[2,1]) + (1/mat[2,2]))
  return(SE)
}