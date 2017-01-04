#' @keywords internal
calcStage2Coefs <- function(Y_i, Zhat_i, V_i, R, q, coefNames){
  
  if(is.null(R)) {
    
    coef_i <- list()
    
    for(i in 1:length(Zhat_i)) { # i <- 10
      coef_i[[i]] <- stats::lm.fit(as.matrix(Zhat_i[[i]]), as.matrix(Y_i[[i]]))$coefficients
    }
    
    coef <- unlist(coef_i)
    names(coef) <- coefNames
    
  } else {
    
    Zhat <- bdiag(Zhat_i)
    
    Y    <- unlist(lapply(Y_i, "["))
    
    W <- rbind2(cbind2(as.matrix(crossprod(Zhat)), t(R)), 
                cbind2(R, matrix(0, nrow(R), nrow(R))))
    
    W <- methods::as(W, "dgeMatrix" )
    
    V <- c(as.numeric(crossprod(Zhat, Y) ), q)
    
    coef <- solve(W, V)[ 1:ncol(Zhat)]
    names(coef) <- coefNames
  }
  
  return(coef)

}