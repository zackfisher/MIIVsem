#' Print method for a MIIV search object
#' 
#' @param x a miivs object
#' @param ... Optional arguments to print, not used by user.
#' 
#' @export
print.miivs <- function(x,...){
  
    z <- x$eqns
  
    for (i in 1:length(z)){
      LHS <- paste(z[[i]]$DVobs, collapse = ", ")
      RHS <- paste(z[[i]]$IVobs, collapse = ", ")
      Instruments <- paste(z[[i]]$MIIVs, collapse = ", ")
      Disturbance <- paste(z[[i]]$CDist, collapse = ", ", sep="")
      modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
      colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
      if (i == 1) {modeqns <- modtemp }
      if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
    }
  
    modeqns$'Composite Disturbance' <- NULL
    
    cat("Model Equation Information \n")
    cat("\n")
    print(
      modeqns,
      quote = FALSE,
      right = FALSE,
      row.names = FALSE,
      print.gap=1
    )
    cat("\n")
    cat("\n")
    

}