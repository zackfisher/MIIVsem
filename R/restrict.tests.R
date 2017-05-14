#'@keywords internal
restrict.tests <- function(x){
  
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
    
    chi.test <- t(R %*% B - q) %*% solve(R %*% coefCov %*% t(R)) %*% (R %*% B - q)
    chi.df   <- nrow(x$r$R)
    chi.p    <- stats::pchisq(chi.test, chi.df, lower.tail = FALSE)
    
    f.test  <- chi.test / nrow(x$r$R)
    f.df1   <- nrow(x$r$R)
    f.df2   <- length(x$eqn)*x$sample.nobs - length(unlist(lapply(x$eqn, "[[","coefficients")))
    f.p     <- stats::pf(f.test, f.df1, f.df2, lower.tail = FALSE)

    return(list(
      chi.test = chi.test, chi.df = chi.df, chi.p = chi.p,
      f.test = f.test, f.df1 = f.df1, f.df2 =f.df2, f.p = f.p)
    )
}