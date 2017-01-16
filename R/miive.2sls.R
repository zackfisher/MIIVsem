#' Two-stage least square estimator for a system of equations
#' 
#' @param d a list containing the system of MIIV estimating equations
#' @param sample.cov sample covariance matrix
#' @param sample.nobs number of observations
#' @param est.only should we only calculate coefficient estimates
#' @param restrictions any equality constraints to be used in estimation

#'@keywords internal
miive.2sls <- function(d, data, sample.cov, sample.mean, sample.nobs, est.only, restrictions){

  # The estimation is done from covariance matrix and mean vector, 
  # so these are calculated first
  
  if(!is.null(data)){
    data <- data[complete.cases(data),]
    sample.cov  <- cov(data)*(nrow(data)-1)/nrow(data)
    sample.nobs <- nrow(data)
    sample.mean <- colMeans(data)
  }
  
  # Build sum of squares and crossproducts matrix (SSCP).
  # From the means, covariances, and n's you can recover the
  # raw sum-of-squares and products matrix for all the variables. 
  # Say the matrix of all the variables is X, with mean vector \bar{x}, 
  # and covariance matrix S, based on sample-size n. Then the SSCP 
  # matrix is X'X = (n - 1)S + n\bar{x}\bar{x}'. You then need to 
  # add the row/column for the constant, which is just n in the 1, 1 
  # position and n\bar{x} elsewhere.

  SSCP <- buildSSCP(sample.cov, sample.nobs, sample.mean)
  
  # Construct the following matrices:
  # XY1: A vector of crossproducts of fitted values from the first stage
  #      regression and the dependent variables. The crossproducts
  #      for all regressions are stacked into one vector
  # XX1: A block diagnonal matrix of of crossproducts of fitted values
  #      from the first stage regressions for all equations.
  # ZV:  A block diagonal matrix of crossproducts between the instuments
  #      and independent variables for all equations
  
  ZV   <- buildBlockDiag(d, SSCP, "IVobs", "MIIVs")
  VV   <- buildBlockDiag(d, SSCP, "MIIVs", "MIIVs", inv = FALSE)
  VY   <- unlist(lapply(d,function(x) SSCP[c("1",x$MIIVs), x$DVobs, drop = FALSE]))

  # DV: Y; EV: Z; MIIVs: V
  # First calculate the part that is used in both equations.
  ZVsVV <- ZV %*% solve(VV)
  XX1   <- ZVsVV %*% t(ZV)
  XY1   <- ZVsVV %*% VY
  
  # TODO: Would it not be more efficient to always calculate the first stage regressions
  # separately than inverting a large VV?
  
  if (is.null(restrictions)){
    
    # b_2sls =  | [Z'V (V'V)^{-1} V'Z]^{-1} |  %*%  |  ZV (V'V)^{-1} V'y  |
    
    coef <- solve(XX1,XY1)
    
  } else {
    
    R <- restrictions$R
    q <- restrictions$q
 
    # b_2sls =  | [Z'V (V'V)^{-1} V'Z]^{-1} | R |       |  ZV (V'V)^{-1} V'y  |
    #           |-------------------------------|  %*%  |---------------------|
    #           |             R             | 0 |       |          q          |
    
    coef <- (solve(rbind(cbind(XX1, t(R)), cbind(R, matrix(0, nrow(R), nrow(R))))) %*%
               rbind(XY1, q))[1:nrow(ZV),]
  }
  
  # TODO: Should the names use Lavaan convetion where regressions of observed 
  # variables on latent variables use =~ instead of x and have the LHS and RHS reversed?
   
  names(coef) <- unlist(lapply(d, function(x) paste0(x$DVlat,"~", c("1", x$IVlat))))
  
  # Start building the return object
  
  res <- list(coefficients = coef,
              sample.cov = sample.cov,
              sample.mean = sample.mean,
              sample.nobs = sample.nobs)
  
  if(!est.only){
    
    # Add coefficients to equations list.
    coefIndex <- unlist(lapply(seq_along(d), function(x) rep(x,(length(d[[x]]$IVobs)+1))))
    coefList  <- split(coef, coefIndex); names(coefList) <- rep("coefficients",length(d))
    d         <- lapply(seq_along(d), function(x) append(d[[x]], coefList[x])) 
    
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
    
    d <- lapply(d, function(eq) { 
      eq$sigma <-  (sample.cov[eq$DVobs, eq$DVobs] +  (t(eq$coefficients[-1]) %*% 
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
    
    res$coefCov <- coefCov

    # Calculate the residual covariance matrix for the full system of
    # equations based on covariance matrix input.
    
    dvs <- unlist(lapply(d, "[[", "DVobs"))
    B   <- diag(length(dvs))
    colnames(B) <- rownames(B) <- dvs
    idx <- do.call("rbind",lapply(d, function(eq) cbind(eq$DVobs, eq$IVobs, eq$coefficients[-1])))
    idx <- idx[idx[,2] %in% dvs, ,drop = FALSE]
    B[idx[,2:1, drop = FALSE]] <- -1*as.numeric(idx[,3])
    
    res$residCov <- t(B) %*% sample.cov[dvs,dvs] %*% B
    
    # Sargan's test from sample covariances (Hayashi, p. 228)
    # TODO: Check for within-equation restrictions 
    #       and alter the df accordingly. What about cross-
    #       equation restrictions, how should this be handled?
    
    # TODO: I passed sample.cov, sample.mean and sample.obs
    #       in the results object. We could move this out 
    #       but then it would be calculated for est.only.
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
    
    res$eqn <- d
    
  }

  return(res)
  
}
