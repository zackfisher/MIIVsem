#' wald test given a model with restrictions
#'@keywords internal
wald <- function(x){
  
    # Greene, 2004 p. 347
  
    # unrestricted coefficients
    d.unr <- est2SLSCoef(
      d = x$eqn, 
      r = NULL, 
      poly.mat = x$sample.polychoric, 
      sscp.mat = x$sample.sscp
    )
    
    d.unr <- est2SLSSigmaSq(
      d = d.unr, 
      cov.mat = x$sample.cov
    )
    
    R <- x$r$R
    q <- x$r$q
    B <- cbind(unlist(lapply(d.unr,"[[","coefficients")))
    
    r.un <- list(R = NULL, q = NULL)

    # unrestricted coefCov
    coefCov <- est2SLSCoefCov( 
      d           = d.unr, 
      poly.mat    = x$sample.polychoric,
      cov.mat     = x$sample.cov,
      mean.vec    = x$sample.mean, 
      acov        = x$asymptotic.cov, 
      acov.sat    = x$asymptotic.cov.sat,
      r           = r.un
    )
    
    W <- t(R %*% B - q) %*% solve(R %*% coefCov %*% t(R)) %*% (R %*% B - q)
    
    wald.df <- nrow(x$r$R)
    wald.p  <- stats::pchisq(W, wald.df, lower.tail = FALSE)

    return(list(wald = W, wald.df = wald.df, wald.p = wald.p))
}