#' build the sum of squares and cross products matrix
#' 
#' Build sum of squares and crossproducts matrix (SSCP).
#' From the means, covariances, and n's you can recover the
#' raw sum-of-squares and products matrix for all the variables. 
#' Say the matrix of all the variables is X, with mean vector bar(x), 
#' and covariance matrix S, based on sample-size n. Then the SSCP 
#' matrix is X'X = (n - 1)S + n bar(x) bar(x)'. You then need to 
#' add the row/column for the constant, which is just n in the 1, 1 
#' position and n bar(x) elsewhere.
#' 
#' @param sample.cov Numeric matrix. A sample variance-covariance matrix. The rownames and colnames must contain the observed variable names.
#' @param sample.mean A sample mean vector.
#' @param sample.nobs Number of observations in the full data frame.
#' 
#'@keywords internal
buildSSCP <- function(sample.cov, sample.mean, sample.nobs){
  
  if (is.null(sample.cov)){
    
    res <- NULL
    
  } else {
    
    res <- matrix(NA, length(sample.mean)+1, length(sample.mean)+1,
                  dimnames = list(c("1", names(sample.mean)),
                                  c("1", names(sample.mean))))
  
    res[1,1]   <- sample.nobs
    res[-1,1]  <- res[1,-1] <- sample.mean*sample.nobs
    res[-1,-1] <- (sample.cov + sample.mean %*% t(sample.mean))*sample.nobs
  }
  
  return(res)
  
}
