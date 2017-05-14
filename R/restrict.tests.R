#' F and Wald test statistics for a set of linear restrictions.
#' 
#' @param x A fitted miive object contraining the model where
#'        at least one restriction is imposted. See the 
#'        \code{\link{miive}} function documentation for
#'        how to impose restrictions using the model 
#'        syntax. 
#' 
#' @details 
#' The internal restrict.tests function provides two test 
#' statistics for a large-sample wald test of any linaer 
#' restrictions imposed on the MIIV-2SLS coefficient matrix.
#' The first statistic is an approximate F and  the second 
#' is Chi-square. Assumptions and  additional details for each
#' test are given by Greene (2003, p. 346-347).
#' 
#' 
#' @return A list containing the following elements:
#'
#' \tabular{ll}{
#' \code{f.test}\tab Wald test statistic\cr
#' \code{f.df}\tab Degrees of freedom for wald test\cr
#' \code{f.p}\tab \cr
#' \code{wald.test}\tab Wald test statistic\cr
#' \code{wald.df}\tab Degrees of freedom for wald test\cr
#' \code{wald.p}\tab P-value associated with Wald test.\cr
#'}
#' 
#' 
#' 
#' @references 
#' 
#' Greene, W. H. (2000). Econometric analysis. Upper Saddle River, N.J: 
#' Prentice Hall.
#' 
#' 
#' @seealso \link{MIIVsem}{miivs}
#' 
#' @keywords MIIV-2SLS MIIV PIV 2sls tsls instrument SEM two-stage least-squares
#'  
#' @export
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
    
    wald.test <- t(R %*% B - q) %*% solve(R %*% coefCov %*% t(R)) %*% (R %*% B - q)
    wald.df   <- nrow(x$r$R)
    wald.p    <- stats::pchisq(wald.test, wald.df, lower.tail = FALSE)
    
    f.test  <- wald.test / nrow(x$r$R)
    f.df1   <- nrow(x$r$R)
    f.df2   <- length(x$eqn)*x$sample.nobs - length(unlist(lapply(x$eqn, "[[","coefficients")))
    f.p     <- stats::pf(f.test, f.df1, f.df2, lower.tail = FALSE)

    return(list(
      wald.test = wald.test, wald.df = wald.df, wald.p = wald.p,
      f.test = f.test, f.df1 = f.df1, f.df2 =f.df2, f.p = f.p)
    )
}