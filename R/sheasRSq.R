#' Calculate Shea's Partial R-squares
#'
#' @param IVobs a vector of endogenous RHS variables.
#' @param MIIVs a vector of instrumental variables.
#' @param sample.cov Numeric matrix. A sample variance-covariance 
#'        matrix. The rownames and colnames attributes must contain all
#'        the observed variable names indicated in the model syntax.
#' @param sample.mean A sample mean vector. If \code{sample.cov} is provided
#'        and the \code{sample.mean} argument is \code{NULL}, intercepts
#'        for all endogenous variables will not be estimated. 
#' @param sample.nobs Number of observations in the full data frame. 
#' 
#'@keywords internal
sheasRSq <- function(IVobs, MIIVs, sample.cov, sample.mean, sample.nobs){
  
  M <- sample.mean
  
  N <- sample.nobs
  
  S <- sample.cov 
  
  sVV <- chol2inv(chol(sample.cov[MIIVs, MIIVs]))

  ZV  <- sample.cov[IVobs , MIIVs, drop = FALSE]

  XX1 <- ZV %*% sVV %*% t(ZV)
  
  if (length(IVobs) == 1){
      
    sheasRSq    <- t(S[IVobs,MIIVs]) %*% solve(S[MIIVs, MIIVs]) %*% 
                     S[IVobs,MIIVs]  %*% solve(S[IVobs, IVobs])
    
    sheasAdjRSq <- 1 - (1 - sheasRSq)*(N-1)/(N - length(MIIVs))
    
    shea <- cbind(sheasRSq, sheasAdjRSq)
    
    rownames(shea) <- IVobs
    
  } else {
    
    shea <- do.call("rbind", lapply(IVobs, function(y){
        
      x2 <- setdiff(IVobs, y)
      
      #M1 <- (XX1[c(y,x2),c(y,x2)] - N * M[c(y,x2)] %*% t(M[c(y,x2)]) ) / (N - 1)
      M1 <- XX1[c(y,x2),c(y,x2)] 
      
      XXt <- S[y,y] -  2*(S[y,x2] %*%  solve(S[x2,x2]) %*% S[x2,y]) +  
             S[y,x2] %*% solve(S[x2,x2]) %*% S[x2,y] 

      XXb <- M1[y,y] - 2*(M1[y,x2] %*% solve(M1[x2,x2]) %*% M1[x2,y]) + 
             M1[y,x2] %*% solve(M1[x2,x2]) %*% M1[x2,y] 
        
      sheasRSq    <- XXb %*% solve(XXt)
      
      sheasAdjRSq <- 1 - (1 - sheasRSq)*(N-1)/(N - length(MIIVs))
        
      sheas <- cbind(sheasRSq, sheasAdjRSq)
      
      rownames(sheas) <- y
      
      sheas
      
    }))
    
  }
  
  colnames(shea) <-c("sheasPartRSq", "sheasPartAdjRSq")
  
  return(shea)
  
}