#' Calculate the minimum eigenvalue statistic of Cragg and Donald (1993)
#'
#' @param IVobs a vector of endogenous RHS variables.
#' @param MIIVs a vector of instrumental variables.
#' @param data matrix or data.frame containing raw data.
#' 
#'@keywords internal
minEigenStat <- function(IVobs, MIIVs, data){
  
  if (is.null(data)){
    stop("Raw data required for minimum eigenvalue statistic.")
  }
  
  # Define operator for taking matrix powers and other functions
  "%^%" <- function(x, n) with(eigen(x), vectors %*% (values^n * t(vectors)))
  M     <- function(A){diag(nrow(A)) - A %*% solve(t(A) %*% A) %*% t(A)}

  if (length(intersect(IVobs, MIIVs)) > 0 ){
    
      X <- as.matrix(cbind(1, data[,intersect(IVobs, MIIVs), drop = FALSE]))
      Z <- as.matrix(data[,MIIVs[!MIIVs %in% intersect(IVobs, MIIVs)], drop = FALSE])
      
  } else {
    
      X <- matrix(1, nrow(data), 1)
      Z <- as.matrix(data[,MIIVs, drop = FALSE])
  }
  
  Y   <- as.matrix(data[,IVobs, drop = FALSE])
  MX1 <- M(X)
  Zb  <- cbind(X, Z)
  MZ  <- M(Zb) 
  Kz  <- ncol(Zb) - 1

  SigVV   <- (t(Y) %*% MZ %*% Y) * (1 / (nrow(data)-Kz))
  SigVV.5 <- SigVV %^% (-0.5)
  G <- t(SigVV.5) %*% t(Y) %*% t(MX1) %*% Z %*% solve(t(Z) %*% MX1 %*% Z) %*% 
       t(Z) %*% MX1 %*% Y %*% (SigVV.5) * (1/Kz)

  minEigen <- min(eigen(G)$values) 
  
  return(minEigen)
}