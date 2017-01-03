#' @keywords internal
calcCoefCovMat <- function(Zhat_i, residCov, R = NULL,q = NULL, numEq, noNA, coefNames){
  
  # residCov <- Omega
  Zhat <- bdiag(Zhat_i)
  class(Zhat) <- "dgCMatrix"
  
  if(class(residCov) != "dspMatrix" ){ class(residCov) <- "dspMatrix" }
  
  sigmaInv <- solve(residCov)
  
  for(i in 1:numEq) {
    for(j in 1:numEq) {
      thisBlock <- sparseMatrix(
        i = which(noNA[noNA[,i],j]),
        j = which(noNA[noNA[,j],i]),
        x = sigmaInv[i,j],
        dims = c(sum(noNA[,i]), sum(noNA[,j])))
      
      if(j == 1) {
        thisRow <- thisBlock
      } else {
        thisRow <- cBind(thisRow, thisBlock)
      }
    }
    
    if( i == 1 ) {
      omegaInv <- thisRow
    } else {
      omegaInv <- rBind(omegaInv, thisRow)
    }
  }
  
  omegaInv <- crossprod(Zhat, omegaInv)
  
  if (is.null(R)) {
    
    coefCov <- solve(omegaInv %*% Zhat)
    
  } else {
    
    W <- rbind2(cbind2(omegaInv %*% Zhat, t(R)),
                cbind2(R, matrix(0, nrow(R), nrow(R))))
    
    coefCov <- as.matrix(solve(W)[1:ncol(Zhat), 1:ncol(Zhat)])
    
  }
  dimnames(coefCov) <- list(coefNames, coefNames)
  return(coefCov)
}