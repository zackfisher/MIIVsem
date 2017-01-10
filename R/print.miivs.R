#' @method print miivs 
#' @export

print.miivs <- function(x,...){
  
  
  for (i in 1:length(x)){
    LHS <- paste(x[[i]]$DVobs, collapse = ", ")
    RHS <- paste(x[[i]]$IVobs, collapse = ", ")
    Instruments <- paste(x[[i]]$MIIVs, collapse = ", ")
    Disturbance <- paste(x[[i]]$CDist, collapse = ", ", sep="")
    modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
    colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
    if (i == 1) {modeqns <- modtemp }
    if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
  }
  
  cat("Model Equation Information \n")
  cat("\n")
  print(modeqns, quote = FALSE, right = FALSE, row.names = FALSE, print.gap=1)
  cat("\n")
  cat("\n")

}