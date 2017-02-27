#' @export
parameterTable <- function(x){
  
    coef.mat <- matrix(x$coefficients, ncol = 1, 
                     dimnames = list(names(x$coefficients),"Estimate"))
  
  if(! is.null(x$coefCov)){
    
    coef.mat <- cbind(coef.mat,
                      Std.Err = sqrt(diag(x$coefCov)),
                      "z-value" = x$coefficients/sqrt(diag(x$coefCov)),
                      "P(>|z|)" = 2*(pnorm(abs(x$coefficients/sqrt(diag(x$coefCov))), lower.tail=FALSE)))
  }
  
  if(! is.null(x$eqn[[1]]$sargan)){
    sarganTests <- do.call(rbind,
                           lapply(x$eqn, function(eq){
                             cbind(c(eq$sargan, rep(NA, length(eq$coefficients)-1)),
                                   c(eq$sargan.df, rep(NA, length(eq$coefficients)-1)),
                                   c(eq$sargan.p, rep(NA, length(eq$coefficients)-1))
                             )
                           }))
    
    colnames(sarganTests) <- c("Sargan", "df", "P(Chi)")
    coef.mat <- cbind(coef.mat, sarganTests)
  }
  
  
  
  if(! is.null(x$varCoefs)){
  
    # Temporary for debugging
    colnames(x$varCoefs) <- c("lhs","rhs","Estimate")
    
    tmp.mat <- matrix(NA, nrow = nrow(varCoefs), ncol = ncol(coef.mat))
    colnames(tmp.mat) <- colnames(coef.mat)
    tmp.mat[,"Estimate"] <- x$varCoefs[,"Estimate"]
    rownames(tmp.at) <- rownames(x$varCoefs)
    coef.mat <- rbind(coef.mat, tmp.mat)
  }
  return(coef.mat)
}