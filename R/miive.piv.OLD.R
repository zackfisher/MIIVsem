#' Polychoric instrumental variable (PIV) estimator for a system of equations
#' 
#' @param d a list containing the system of MIIV estimating equations
#' @param sample.cov sample covariance matrix
#' @param sample.nobs number of observations
#' @param est.only should we only calculate coefficient estimates
#' @param restrictions any equality constraints to be used in estimation
#' @param piv.opts Options to pass to lavaan's \code{\link[lavCor]{lavCor}} function.

#'@keywords internal
miive.piv <- function(d, data, sample.cov, sample.mean, sample.nobs, est.only, 
                      restrictions, piv.opts = piv.opts){

  # The regression coefficients are estimated from
  # using the polychoric correlation matrix, so this must 
  # be constructed first. The raw data is needed so throw
  # error here if not available. 
  if(is.null(data)){
    stop(paste(
      "miive: Raw data is required to generate polychoric correlation matrix.")
    )
  }
  
  if(!is.null(data)){
    
    # NOTE: We could use missing = "FIML" for the PIV estimator when
    # generating the polychoric correlation matrix and acov.
    # data <- data[complete.cases(data),]
    
    # TODO: Need to identify any exogenous factors. These must be
    # handled seperately. How does lavaan handle this in generation
    # of polychoric correlation matrix. 
    
    pcr  <- unclass(lavaan::lavCor(
      data, 
      output= "cor", 
      #se = piv.opts["se"], 
      estimator = piv.opts["estimator"],
      ordered = colnames(data)[!apply(data,2,is.numeric)]
    ))
    
    
    # Generate the asymptotic covariance matrix of polychoric 
    # correlations.
    
    acov <- unclass(lavaan::vcov(lavaan::lavCor(
      data, 
      output = "fit", 
      se = piv.opts["se"] , 
      estimator = piv.opts["estimator"],
      ordered = colnames(data)[!apply(data,2,is.numeric)]
    ))) 
    
    # Remove thresholds from acov.pcr
    acov <- acov[1:(1/2*nrow(pcr)*(nrow(pcr)-1)),1:(1/2*nrow(pcr)*(nrow(pcr)-1))]
    
    sample.nobs <- nrow(data)
    sample.mean <- NULL
  }
  
  # Construct the following matrices:
  # XY1: A vector of crossproducts of fitted values from the first stage
  #      regression and the dependent variables. The crossproducts
  #      for all regressions are stacked into one vector
  # XX1: A block diagnonal matrix of of crossproducts of fitted values
  #      from the first stage regressions for all equations.
  # ZV:  A block diagonal matrix of crossproducts between the instuments
  #      and independent variables for all equations
  
  ZV   <- buildBlockDiag(d, pcr, "IVobs", "MIIVs", pcr = TRUE)
  VV   <- buildBlockDiag(d, pcr, "MIIVs", "MIIVs", pcr = TRUE)
  VY   <- unlist(lapply(d,function(x) pcr[x$MIIVs, x$DVobs, drop = FALSE]))

  # DV: Y; EV: Z; MIIVs: V
  # First calculate the part that is used in both equations.
  ZVsVV <- ZV %*% chol2inv(chol(VV))
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
 
    #b_2sls =  | [Z'V (V'V)^{-1} V'Z]^{-1} | R |       |  ZV (V'V)^{-1} V'y  |
    #          |-------------------------------|  %*%  |---------------------|
    #          |             R             | 0 |       |          q          |
    # as.numeric makes coef a vector instead of a matrix
    coef <- as.numeric((solve(rbind(cbind(XX1, t(R)), cbind(R, matrix(0, nrow(R), nrow(R))))) %*%
              rbind(XY1, q))[1:nrow(ZV),])
  }
  
  # TODO: Should the names use Lavaan convetion where regressions of observed 
  # variables on latent variables use =~ instead of x and have the LHS and RHS reversed?
   
  names(coef) <- unlist(lapply(d, function(x) paste0(x$DVlat,"~", x$IVlat)))
  
  # Start building the return object
  
  res <- list(coefficients = coef,
              sample.cov   = pcr,
              sample.mean  = sample.mean,
              sample.nobs  = sample.nobs)

  # Add coefficients to equations list.
  coefIndex <- unlist(lapply(seq_along(d), function(x) rep(x,(length(d[[x]]$IVobs)))))
  coefList  <- split(coef, coefIndex); names(coefList) <- rep("coefficients",length(d))
  d         <- lapply(seq_along(d), function(x) append(d[[x]], coefList[x])) 
  
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
    
    # d <- lapply(d, function(eq) { 
    #   eq$sigma <-  (sample.cov[eq$DVobs, eq$DVobs] +  (t(eq$coefficients[-1]) %*% 
    #                sample.cov[c(eq$IVobs), c(eq$IVobs)] %*% eq$coefficients[-1]) -
    #                (2 * sample.cov[eq$DVobs, c(eq$IVobs)] %*% eq$coefficients[-1])) 
    #   eq
    # })
    
    #sig <- diag(unlist(lapply(d, function(eq) rep(eq$sigma, length(eq$coefficients)))))
    
    if (is.null(restrictions)){
      # Standard Errors for the PIV Regression Coefficients
      # For each equation in d:
      #   Obtain partial derivatives of \theta with respect to the polychoric
      #   correlations among the instrumental variables, the polychoric 
      #   correlations between the instruments and explanatory variables,
      #   and between the instruments and the dependent variable. 
      #
      # With those derivatives in hand we must assemble the K matrix to 
      # correspond in structure to the asymptotic covariance matrix of 
      # the polychoric correlations. Eq. 25, page 315. Bollen, K. A., 
      # & Maydeu-Olivares, A. (2007). A Polychoric Instrumental Variable 
      # (PIV) Estimator for Structural Equation Models with Categorical 
      # Variables. Psychometrika, 72(3), 309–326. 

      K       <- buildKmatrix(d, pcr)
      coefCov <- K %*% acov %*% t(K) 
      
      
      
    } else {
      
      # TODO: It is likely restruction need to be handled in the
      #       calculation of the partial derivatives. 
      
      
      #R0 <- matrix(0, ncol=nrow(R), nrow=nrow(R))
      #coefCov <- solve(rbind(cbind((XX1 %*% t(solve(sig))), t(R)), 
      #                       cbind(R, R0)))[1:nrow(XX1), 1:nrow(XX1)]
    }
    
    res$coefCov <- coefCov

    # Calculate the residual covariance matrix for the full system of
    # equations based on covariance matrix input.
    
    dvs <- unlist(lapply(d, "[[", "DVobs"))
    B   <- diag(length(dvs))
    colnames(B) <- rownames(B) <- dvs
    idx <- do.call("rbind",lapply(d, function(eq) cbind(eq$DVobs, eq$IVobs, eq$coefficients)))
    idx <- idx[idx[,2] %in% dvs, ,drop = FALSE]
    B[idx[,2:1, drop = FALSE]] <- -1*as.numeric(idx[,3])
    
    res$residCov <- t(B) %*% pcr[dvs,dvs] %*% B
    
    # Sargan's test from sample covariances (Hayashi, p. 228)
    # TODO: Check for within-equation restrictions 
    #       and alter the df accordingly. What about cross-
    #       equation restrictions, how should this be handled?
    
    # TODO: I passed sample.cov, sample.mean and sample.obs
    #       in the results object. We could move this out 
    #       but then it would be calculated for est.only.
    
    # TODO: Is this valid when using the polychoric correlations?
    #
    # d <- lapply(d, function(eq) { 
    #   eq$sargan <-  (
    #     t(sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] - 
    #         sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*% 
    #         eq$coefficients[-1]) %*% 
    #       solve(sample.cov[eq$MIIVs,eq$MIIVs]) %*% 
    #       (sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] - 
    #          sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*% 
    #          eq$coefficients[-1]) /  eq$sigma)*sample.nobs
    #   eq$sargan.df <- length(eq$MIIVs) - length(eq$IVobs)
    #   eq$sargan.p <- pchisq(eq$sargan, eq$sargan.df, lower.tail = FALSE)
    #   eq
    # })
  }

  res$eqn <- d
  
  return(res)
  
}
