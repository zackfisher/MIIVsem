#' Two-stage least square estimator for a system of equations
#' 
#' @param d a list containing the system of MIIV estimating equations
#' @param sample.cov sample covariance matrix
#' @param sample.nobs number of observations
#' @param se should variance covariance matrix of the estimates be calculated
#' @param restrictions any equality constraints to be used in estimation

#'@keywords internal
miive.2sls <- function(d, data, sample.cov, sample.mean, sample.nobs, se, restrictions){

  # The estimation is done from covariance matrix and mean vector, so these are calculated first
  
  if(!is.null(data)){
    data <- data[complete.cases(data),]
    sample.cov <- cov(data)*(nrow(data)-1)/nrow(data)
    sample.nobs <- nrow(data)
    sample.mean <- colMeans(data)
  }
  
  
  # TODO: Explain what SSP, ZV, VV, and VY are
  
  SSP <- buildSSP(sample.cov, sample.nobs, sample.mean)
  ZV  <- buildBlockDiag(d, SSP, "IVobs", "MIIVs")
  VV  <- buildBlockDiag(d, SSP, "MIIVs", "MIIVs")
  VY  <- unlist(lapply(d,function(x) SSP[c("1",x$MIIVs), x$DVobs, drop = FALSE]))

  # DV: Y; EV: Z; MIIVs: V
  XX1 <- ZV %*% solve(VV) %*% t(ZV)
  XY1 <- ZV %*% solve(VV) %*% VY
  
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
  
  # TODO: Should the names use Lavaan convention where regressions of observed 
  # variables on latent variables use =~ instead of x and have the LHS and RHS reversed? 
  
  rownames(coef) <- unlist(lapply(d, function(x) paste0(x$DVlat,"~", c("1", x$IVlat))))
  
  # Add coefficients to equations list.
  coefIndex <- unlist(lapply(seq_along(d), function(x) rep(x,(length(d[[x]]$IVobs)+1))))
  coefList  <- split(coef, coefIndex); names(coefList) <- rep("coefficients",length(d))
  d         <- lapply(seq_along(d), function(x) append(d[[x]], coefList[x])) 
  
  # Start building the return object
  
  res <- list(coefficients = coef)
  
  if(se){
    
    #         | sigma_{11}                   |
    # Sigma = |            ...               |
    #         |                  sigma_{neq} | 
    #
    # sigma_{11} = S_{YY} + b1' S_{ZZ} b1 - 2 * S_{YZ} b1
    # b1 does not contain intercepts.
    
    Sigma <- lavaan::lav_matrix_bdiag(lapply(d, function(eq){
      (sample.cov[eq$DVobs, eq$DVobs] +  (t(eq$coefficients[-1]) %*% 
      sample.cov[c(eq$IVobs), c(eq$IVobs)] %*% eq$coefficients[-1]) -
         (2 * sample.cov[eq$DVobs, c(eq$IVobs)] %*% eq$coefficients[-1])) 
    }))
    
    sig <- diag(unlist(lapply(seq_along(d), function(i) rep(Sigma[i,i], length(d[[i]]$coefficients)))))
    
    res$Sigma <- Sigma
    
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
    diag(B) <- 1
    
    res$ResidCov <- t(B) %*% sample.cov[dvs,dvs] %*% B

    res$eqn <- d
    
  }

  return(res)
  
}
