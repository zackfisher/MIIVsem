#' @method print miive 
#' @export
print.miive <- function(x,  digits = max(3, getOption("digits") - 2),...){
  
  header.mat <- as.data.frame(rbind(
    c("Number of observations: ", x$sample.nobs),
    c("Number of equations: ", length(x$eqn)),
    c("Estimator: ", x$estimator)), row.names = NULL)
  
  for (i in 1:ncol(header.mat)) {
    header.mat[, i] = format(header.mat[, i], 
                            justify = ifelse(i == 1, "left", "right"))
  }
  
  cat("\n")
  cat(paste0("MIIVsem (", packageVersion("MIIVsem"),") results"), "\n")
  cat("\n")
  print(header.mat,quote = FALSE,row.names = FALSE, col.names = FALSE)
  #cat("Number of observations:", x$sample.nobs, "\n")
  #cat("Number of equations:", length(x$eqn), "\n")
  #cat("Estimator:", x$estimator, "\n")
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