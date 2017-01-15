# Miscellanneous helper functions
# A helper function for binding a list of matrices by their names
#@param x a list of matrices
#@param empty value to be used for empty cells
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


buildBlockDiag <- function(d, mat, row, col, inv = FALSE){
 
  res <- lavaan::lav_matrix_bdiag(
    if (inv){
      lapply(d, function(eq){
        solve(mat[c("1",eq[[row]]), c("1",eq[[col]]), drop = FALSE]) 
      })
    } else {
      lapply(d, function(eq){
        mat[c("1",eq[[row]]), c("1",eq[[col]]), drop = FALSE] 
      })
    })
}

buildSSCP <- function(sample.cov, sample.nobs, sample.means){
  
  res <- matrix(NA, length(sample.means)+1, length(sample.means)+1,
                dimnames = list(c("1", names(sample.means)),
                                   c("1", names(sample.means))))
  
  res[1,1]   <- sample.nobs
  res[-1,1]  <- res[1,-1] <- sample.means*sample.nobs
  res[-1,-1] <- (sample.cov + sample.means %*% t(sample.means))*sample.nobs
  return(res)
}


