#' Calculate Shea's Partial R-squares
#'
#' @param object a fitted miive object.
#' 
#'@keywords internal
sheasRSq <- function(IVobs, MIIVs, XX1, sample.cov, sample.mean, sample.nobs){
  
  M <- sample.mean
  N <- sample.nobs
  # undo smaple.cov rescale to match Stata results
  S <- sample.cov * N/(N-1)
  
  if (length(IVobs) == 1){
      
    sheasRSq    <- t(S[IVobs,MIIVs]) %*% solve(S[MIIVs, MIIVs]) %*% 
                     S[IVobs,MIIVs] %*% solve(S[IVobs, IVobs])
    
    sheasAdjRSq <- 1 - (1 - sheasRSq)*(N-1)/(N - length(MIIVs))
    shea <- cbind(sheasRSq, sheasAdjRSq)
    rownames(shea) <- y
    
  } else {
    
    shea <- do.call("rbind", lapply(IVobs, function(y){
        
      x2 <- setdiff(IVobs, y)
      
      M1 <- (XX1[c(y,x2),c(y,x2)] - N * M[c(y,x2)] %*% t(M[c(y,x2)]) ) / (N - 1)
      
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