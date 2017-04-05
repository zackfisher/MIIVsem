#' @method print miive 
#' @export
print.miive <- function(x,...){

  # MIIVsem version number
  cat(paste0("MIIVsem (", packageVersion("MIIVsem"),") results"), "\n\n")
  
  w1 <- 30 # width of column 1
  w2 <- 30 # width of column 2
  
  head.txt  <- do.call("rbind",
               list(c("Number of observations", x$sample.nobs),
                    c("Number of equations", length(x$eqn)),
                    c("Estimator", x$estimator)))
  for(i in 1:nrow(head.txt)){
    cat(sprintf("%-*s %*s\n", w1, head.txt[i,1], w2, head.txt[i, 2]));
  }; cat("\n")
  
  # This is only temporary for debugging.
  # 
  # TODO: switch to lavaan-style encoding, e.g. 
  #       RHS =~ LHS for latent variable and similar
  #       output when summary command is called.
  
  coef.mat <- parameterTable(x)
  
  coef.mat.df <- as.data.frame(coef.mat)
  coef.mat <- as.matrix(format.data.frame(coef.mat.df, digits = 3, na.encode = FALSE))
  coef.mat[is.na(coef.mat.df)] <- NA
  
  print(coef.mat, digits = 3, na.print = "", quote = FALSE, justify = "none")
  
  if(! is.null(x$varCoefs)){
    
    cat("\n")
    cat(paste0("Variance and Covariance Parameter Estimates"), "\n")
    cat("\n")
    # Temporary for debugging
    colnames(x$varCoefs) <- c("lhs","rhs","Estimate")
    print(x$varCoefs[,3,drop=FALSE], digits = digits, na.print = "", quote = FALSE, justify = "none")
  }
}