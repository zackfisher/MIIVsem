#' Calculate Shea's Partial R-squares
#'
#' @param object a fitted miive object.
#' 
#'@keywords internal
sheasRSq <- function(object){
  
  M <- object$sample.mean
  N <- object$sample.nobs
  # undo smaple.cov rescale to match Stata results
  S <- object$sample.cov * N/(N-1)
  
  object$eqn <- lapply(x$eqn, function(eq){
    
    if (length(eq$IVobs) == 1){
      
        sheasRSq    <- t(S[eq$IVobs,eq$MIIVs]) %*% solve(S[eq$MIIVs, eq$MIIVs]) %*% 
                         S[eq$IVobs,eq$MIIVs] %*% solve(S[eq$IVobs, eq$IVobs])
        sheasAdjRSq <- 1 - (1 - sheasRSq)*(N-1)/(N - length(eq$MIIVs))
        shea <- cbind(sheasRSq, sheasAdjRSq)
        rownames(shea) <- y
        shea
    
    } else {
    
      sheas <- do.call("rbind", lapply(eq$IVobs, function(y){
        
        x2 <- setdiff(eq$IVobs, y)
      
        M1 <- (eq$XX1[c(y,x2),c(y,x2)] - N * M[c(y,x2)] %*% t(M[c(y,x2)]) ) / (N - 1)
      
        XXt <- S[y,y] -  2*(S[y,x2] %*%  solve(S[x2,x2]) %*% S[x2,y]) +  
               S[y,x2] %*% solve(S[x2,x2]) %*% S[x2,y] 

        XXb <- M1[y,y] - 2*(M1[y,x2] %*% solve(M1[x2,x2]) %*% M1[x2,y]) + 
               M1[y,x2] %*% solve(M1[x2,x2]) %*% M1[x2,y] 
        
        sheasRSq    <- XXb %*% solve(XXt)
        sheasAdjRSq <- 1 - (1 - sheasRSq)*(N-1)/(N - length(eq$MIIVs))
        shea <- cbind(sheasRSq, sheasAdjRSq)
        rownames(shea) <- y
        shea
      }))
      
      colnames(sheas) <-c("sheasPartRSq", "sheasPartAdjRSq")
    
    }
  
  })
  
  return(object)
}