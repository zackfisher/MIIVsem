# Miscellanneous helper functions

# A helper function for binding a list of matrices by their names
#'@param x a list of matrices
#'@param empty value to be used for empty cells
#'@keywords internal

miivhelper.bindmatrices <- function(x, empty){
  
  coln <- unique(unlist(lapply(x,colnames)))
  rown <- unique(unlist(lapply(x,rownames)))
  
  mat <- matrix(empty,length(rown), length(coln),
                dimnames = list(rown,coln))
  
  for(i in seq_along(x)){
    mat[rownames(x[[i]]), colnames(x[[i]])] <- x[[i]]
  }
  
  mat
}
