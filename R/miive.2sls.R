#' two-stage least square estimator for a system of equations
#' 
#' @param d a list containing equation information
#' @param g a list containing data matrices and characteristics
#' @param r a list containing coefficient restrictions 
#' @param est.only should we only calculate coefficient estimates
#' @param se se estimation
#' 
#'@keywords internal
miive.2sls <- function(d, g, r, est.only, se, missing, var.cov, sarg.adjust="none"){
  
  #------------------------------------------------------------------------#
  # MIIV-2SLS point estimates
  #------------------------------------------------------------------------#
  d <- est2SLSCoef(d, r, g$sample.polychoric, g$sample.sscp)
  coefficients <- unlist(lapply(d,"[[","coefficients"))
  names(coefficients) <- unlist(lapply(d,function(eq)names(eq$coefficients)))
  
  coefCov <- NULL
  
  #------------------------------------------------------------------------#
  # MIIV-2SLS se 
  #------------------------------------------------------------------------#
  if(!est.only){ 
    
    coef.names <- unlist(lapply(d, function(eq) { 
      paste0(eq$DVlat, "~",if(eq$categorical) eq$IVlat else c("1",eq$IVlat))
    }))
    
    # calculate sigma^2
    d <- est2SLSSigmaSq(d, cov.mat = g$sample.cov)
  
     if(se != "boot" & se != "bootstrap" & missing != "twostage"){ 
        coefCov <- est2SLSCoefCov( 
          d           = d, 
          poly.mat    = g$sample.polychoric,
          cov.mat     = g$sample.cov,
          mean.vec    = g$sample.mean, 
          acov        = g$asymptotic.cov, 
          acov.sat    = g$asymptotic.cov.sat,
          r           = r 
        )
        
     }
  
  #------------------------------------------------------------------------#
  # MIIV-2SLS overidentification tests
  #------------------------------------------------------------------------#

    # Sargan's test (Hayashi, p. 228)
    d <- lapply(d, function(eq) {
      
      #[ t(Svy - Svz * B) %*% SvvInv %*% (Svy - Svz * B) / sig ] * N
      
      if (eq$categorical | (length(eq$MIIVs) - length(eq$IVobs)) < 1 ){

        eq$sargan <- NA; eq$sargan.df <- NA; eq$sargan.p  <- NA

      } else {

        eq$sargan.df <- length(eq$MIIVs) - length(eq$IVobs)

        eq$sargan <-
          (
            t( g$sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] -
               g$sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*%
               eq$coefficients[-1]) %*%
               solve(g$sample.cov[eq$MIIVs,eq$MIIVs] )
            %*%
             ( g$sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] -
               g$sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*%
               eq$coefficients[-1] )
            /  eq$sigma
          ) * g$sample.nobs

          eq$sargan.p <- stats::pchisq(

            eq$sargan, eq$sargan.df, lower.tail = FALSE

          )
      }
      eq
     })
    
    if (sarg.adjust != "none"){
      sarg.pvals <- unlist(lapply(d, "[[", "sargan.p"))
      adj.sarg.pvals <- stats::p.adjust(sarg.pvals, sarg.adjust)
      
      for (i in 1:length(adj.sarg.pvals)){
        d[[i]]$sargan.p <- adj.sarg.pvals[i]
      }
      
    }
    
  }
  
  g$coefficients <- coefficients
  g$coefCov      <- coefCov
  g$eqn          <- d
  
  return(g)
  
}

