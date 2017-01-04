#' @keywords internal
calcResidCovMat <- function(resids_i, numEq, noNA, oneSigma = FALSE,
                            centered = FALSE, diagOnly = TRUE) {
  
  if (!oneSigma){
    if (any(is.na(noNA))) {
      resids <- matrix(NA, nrow(noNA), ncol(noNA))
      for (i in 1:numEq) {resids[noNA[,i],i] <- resids_i}
    } else {
    }
    
    if(centered) {
      for(i in 1:numEq) {
        resids[,i] <- resids[,i] - mean(resids[noNA[,i],i])
      }
    }
    
    Omega <- matrix(0, numEq, numEq)
    
    noNARow <- rowSums(!noNA) == 0
    sumRows <- sum(noNARow)
    
    if (diagOnly){
      for(i in 1:numEq) {
        Omega[i, i] <- sum(resids[noNARow, i]^2) / sumRows 
      }
    }
    
    if (!diagOnly){  
      for(i in 1:numEq) {
        for(j in 1:numEq) {
          Omega[i,j] <-sum(resids[noNARow, i] * resids[noNARow, j]) / sumRows 
        }
      }
      Omega <- kronecker(Omega, diag())
    }
  }
  
#  Mikko: I commented this out because the Y_i is no defined. This leads to a warning
#  during build and would produce and error if the code is executed.
#
#  if(oneSigma){
#    resids <- do.call("rbind", resids_i)
#    Omega  <- sum( resids^2 ) / length(unlist(Y_i))
#    Omega  <- Omega * diag(length(unlist(Y_i)))
#  }
  
  Omega <- methods::as(Omega, "dspMatrix")
  return(Omega)
}