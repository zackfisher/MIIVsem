#' Two-stage least square estimator for a system of equations
#' 
#' @param d a list containing the system of MIIV estimating equations
#' @param sample.cov sample covariance matrix
#' @param sample.nobs number of observations
#' @param est.only should we only calculate coefficient estimates
#' @param restrictions any equality constraints to be used in estimation
#' @param piv.opts Options to pass to lavaan's \code{\link[lavCor]{lavCor}} function.
#' @param factorIndex logical vector indexing whether variables are factors

#'@keywords internal
miive.2sls <- function(d, data, sample.cov, sample.mean, sample.nobs, est.only, piv.opts, factorIndex){

  # Build sum of squares and crossproducts matrix (SSCP).
  # From the means, covariances, and n's you can recover the
  # raw sum-of-squares and products matrix for all the variables. 
  # Say the matrix of all the variables is X, with mean vector \bar{x}, 
  # and covariance matrix S, based on sample-size n. Then the SSCP 
  # matrix is X'X = (n - 1)S + n\bar{x}\bar{x}'. You then need to 
  # add the row/column for the constant, which is just n in the 1, 1 
  # position and n\bar{x} elsewhere.
  
  if(any(factorIndex)){
    
    pcr    <- TRUE
    
    estMat <- unclass(lavaan::lavCor(data, output= "cor", 
                      estimator = piv.opts["estimator"],
                      ordered = colnames(data)[factorIndex]))
    
    # Generate the asymptotic covariance matrix of polychoric 		
    # correlations.		
    acov <- unclass(lavaan::vcov(lavaan::lavCor(data, output = "fit", 		
            se = piv.opts["se"] ,estimator = piv.opts["estimator"],		
            ordered = colnames(data)[factorIndex]	))) 
    
    # Remove thresholds from acov.pcr		
    acov <- acov[1:(1/2*nrow(estMat)*(nrow(estMat)-1)),1:(1/2*nrow(estMat)*(nrow(estMat)-1))]
    
    
  
  } else {
    
    pcr    <- FALSE
    
    estMat <- buildSSCP(sample.cov, sample.nobs, sample.mean)
    
  }
  
  #-------------------------------------------------------#  
  # Build Restriction Matrices.
  # returns NULL if there are no restrictions,
  # otherwise returns a list containing the 'R' matrix
  # and 'q' vector, as well as a vector 'cons' of 
  # the constrained coefficients.
  #
  # We throw an error here if there is a restriction 
  # between a equation with d$categorical = TRUE and
  # d$categorical = FALSE.
  #-------------------------------------------------------#  
  restrictions <- buildRestrictMat(d, pcr)
  #-------------------------------------------------------#  

  # Construct the following matrices:
  # XY1: A vector of crossproducts of fitted values from the first stage
  #      regression and the dependent variables. The crossproducts
  #      for all regressions are stacked into one vector
  # XX1: A block diagnonal matrix of of crossproducts of fitted values
  #      from the first stage regressions for all equations.
  # ZV:  A block diagonal matrix of crossproducts between the instuments
  #      and independent variables for all equations
  ZV   <- buildBlockDiag(d, estMat, "IVobs", "MIIVs", inv = FALSE, pcr)
  VV   <- buildBlockDiag(d, estMat, "MIIVs", "MIIVs", inv = FALSE, pcr)
  VY   <- unlist(lapply(d, function(x){estMat[if(pcr) x$MIIVs else c("1",x$MIIVs), x$DVobs, drop = FALSE]}))
    

  # DV: Y; EV: Z; MIIVs: V
  # First calculate the part that is used in both equations.
  #ZVsVV <- ZV %*% chol2inv(chol(VV))
  ZVsVV <- lavaan::lav_matrix_bdiag(lapply(d, function(eq){
    estMat[c("1",eq[["IVobs"]]), c("1",eq[["MIIVs"]]), drop = FALSE] %*%
      chol2inv(chol(estMat[c("1",eq[["MIIVs"]]), c("1",eq[["MIIVs"]]), drop = FALSE]))
  }))
  
  XX1   <- ZVsVV %*% t(ZV)
  XY1   <- ZVsVV %*% VY
  
  # TODO: Would it not be more efficient to always calculate the first stage regressions
  # separately than inverting a large VV? This should be checked with large models
  # NOTE: In the small bit of profiling I did with splitting this up by
  # equation it didn't make a big difference with the political democracy
  # example.  It might though when there are more observed variables in the model,
  # leading to more instruments per equation. Maybe toggling between the two
  # depending on the number of questions. 
  
  if (is.null(restrictions)){
    
    # b_2sls =  | [Z'V (V'V)^{-1} V'Z]^{-1} |  %*%  |  ZV (V'V)^{-1} V'y  |
    # as.numeric makes coef a vector instead of a matrix
    coef <- as.numeric(solve(XX1,XY1))
    
  } else {
    
    R <- restrictions$R
    q <- restrictions$q
 
    # b_2sls =  | [Z'V (V'V)^{-1} V'Z]^{-1} | R |       |  ZV (V'V)^{-1} V'y  |
    #           |-------------------------------|  %*%  |---------------------|
    #           |             R             | 0 |       |          q          |
    # as.numeric makes coef a vector instead of a matrix
    coef <- as.numeric((solve(rbind(cbind(XX1, t(R)), cbind(R, matrix(0, nrow(R), nrow(R))))) %*%
               rbind(XY1, q))[1:nrow(ZV),])
  }
  
  # TODO: Should the names use Lavaan convetion where regressions of observed 
  # variables on latent variables use =~ instead of x and have the LHS and RHS reversed?
   
  names(coef) <- unlist(lapply(d, function(x) {
    paste0(x$DVlat,"~", if(pcr) x$IVlat else c("1",x$IVlat) )
  }))
  
  # Start building the return object
  
  res <- list(coefficients = coef,
              sample.cov = sample.cov,
              sample.mean = sample.mean,
              sample.nobs = sample.nobs,
              restrictions = restrictions,
              estimator = ifelse(pcr, "MIIV-2SLS (PIV)", "MIIV-2SLS"))

  # Add coefficients to equations list.
  coefIndex <- unlist(lapply(seq_along(d), function(x) {
    rep(x,(length(d[[x]]$IVobs) + ifelse(pcr, 0, 1)))
  }))
  
  coefList  <- split(coef, coefIndex); 
  names(coefList) <- rep("coefficients",length(d))
  d  <- lapply(seq_along(d), function(x) append(d[[x]], coefList[x])) 
  
  if(!est.only){
    
    
    #         | sigma_{11}                   |
    # Sigma = |            ...               |
    #         |                  sigma_{neq} | 
    #
    # sigma_{11} = S_{YY} + b1' S_{ZZ} b1 - 2 * S_{YZ} b1
    # b1 does not contain intercepts.
    
    # Add equation-level sigma^2 to eqns list instead of block diagonal 
    # big Sigma to results which we obtain in full later. This makes 
    # calculation of equation level statistics easier. 
    
    # TODO: How should equation-level sigma^2 be handled when cross-
    #       equation restrictions are present?
    
    if (pcr){ # begin PCR

      K       <- buildKmatrix(d, estMat)
      coefCov <- K %*% acov %*% t(K) 


    } else { # begin SSCP
      
      d <- lapply(d, function(eq) { 
        eq$sigma <-(sample.cov[eq$DVobs, eq$DVobs] + (t(eq$coefficients[-1]) %*% 
                    sample.cov[c(eq$IVobs), c(eq$IVobs)] %*% eq$coefficients[-1]) -
                   (2 * sample.cov[eq$DVobs, c(eq$IVobs)] %*% eq$coefficients[-1])) 
        eq
      })
      
      sig <- diag(unlist(lapply(d, function(eq) rep(eq$sigma, length(eq$coefficients)))))
      
      if (is.null(restrictions)){
        
        coefCov <- solve(XX1 %*% t(solve(sig)))
        
      } else {
        
        R0 <- matrix(0, ncol=nrow(R), nrow=nrow(R))
        coefCov <- solve(rbind(cbind((XX1 %*% t(solve(sig))), t(R)), 
                               cbind(R, R0)))[1:nrow(XX1), 1:nrow(XX1)]
      }
    }
    
    
    res$coefCov <- coefCov

    # Calculate the residual covariance matrix for the full system of
    # equations based on covariance matrix input.
    
    dvs <- unlist(lapply(d, "[[", "DVobs"))
    B   <- diag(length(dvs))
    colnames(B) <- rownames(B) <- dvs
    idx <- do.call("rbind",lapply(d, function(eq){
      cbind(eq$DVobs, eq$IVobs, if (pcr) eq$coefficients else eq$coefficients[-1] )
    }))
    idx <- idx[idx[,2] %in% dvs, ,drop = FALSE]
    B[idx[,2:1, drop = FALSE]] <- -1*as.numeric(idx[,3])
    
    res$residCov <- t(B) %*% (if (pcr) estMat else sample.cov)[dvs,dvs] %*% B
    
    # Sargan's test from sample covariances (Hayashi, p. 228)
    # TODO: Check for within-equation restrictions 
    #       and alter the df accordingly. What about cross-
    #       equation restrictions, how should this be handled?
    
    # TODO: I passed sample.cov, sample.mean and sample.obs
    #       in the results object. We could move this out 
    #       but then it would be calculated for est.only.
    
    if (pcr){
      # TODO: Theoretical justification needed for Sargan based
      #       on polychoric correlations.
      d <- lapply(d, function(eq) { 
        eq$sargan    <- NULL
        eq$sargan.df <- NULL
        eq$sargan.p  <- NULL
        eq
      })
      
    } else {
    d <- lapply(d, function(eq) { 
      eq$sargan <-  (
        t(sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] - 
            sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*% 
            eq$coefficients[-1]) %*% 
          solve(sample.cov[eq$MIIVs,eq$MIIVs]) %*% 
          (sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] - 
             sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*% 
             eq$coefficients[-1]) /  eq$sigma)*sample.nobs
      eq$sargan.df <- length(eq$MIIVs) - length(eq$IVobs)
      eq$sargan.p <- pchisq(eq$sargan, eq$sargan.df, lower.tail = FALSE)
      eq
     })
   }

  res$eqn <- d
  return(res)
  }
}
