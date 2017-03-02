#' @method print miive 
#' @export
print.miive <- function(x,  digits = max(3, getOption("digits") - 2),...){
  
  cat("\n")
  cat(paste0("MIIVsem (", packageVersion("MIIVsem"),") results"), "\n")
  cat("\n")
  cat("Number of observations:", x$sample.nobs, "\n")
  cat("Number of equations:", length(x$eqn), "\n")
  cat("\n")
  
  # This is only temporary for debugging.
  # 
  # TODO: switch to lavaan-style encoding, e.g. 
  #       RHS =~ LHS for latent variable and similar
  #       output when summary command is called.
  
  coef.mat <- parameterTable(x)
  
  print(coef.mat, digits = digits, na.print = "", quote = FALSE, justify = "none")
  
  if(! is.null(x$varCoefs)){
    
    cat("\n")
    cat(paste0("Variance and Covariance Parameter Estimates"), "\n")
    cat("\n")
    # Temporary for debugging
    colnames(x$varCoefs) <- c("lhs","rhs","Estimate")
    print(x$varCoefs[,3,drop=FALSE], digits = digits, na.print = "", quote = FALSE, justify = "none")
  }
}