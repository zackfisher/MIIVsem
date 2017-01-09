#' Two-stage least square estimator
#'
#' The estimator handles missing data by equationwise deletion, which 
#' means that all equations are estimated using cases that have all data
#' for that equation. It is possible that different equations 
#' are estimated with different subsamples. If other kind of 
#' missing data prcessing is needed (e.g.) listwise deletion,
#' that must be done prior to calling this function
#' 
#' @param d list of MIIV equations
#' @param data optional data.frame
#' @param sample.cov optional sample covariance matrix
#' @param se should variance covariance matrix of the estimates be calculated
#' 
#'@keywords internal

miive.2sls <- function(d, data, sample.cov, sample.mean, sample.nobs, se, restrictions){
  
  # Estimate the system
  
  # If we have data, then calculate a sample covariance matrix for this equation
  if(!is.null(data)){

      data <- data[complete.cases(data),]
      
      results <- miive.2sls.system(d, cov(data)*(nrow(data)-1)/nrow(data), colMeans(data), nrow(data), se, restrictions)
  }
    
  # Else use the sample covariance matrix given as argument
  else{
    
      results <- miive.2sls.system(d, sample.cov, sample.mean, sample.nobs, se, restrictions)
      
  }    

  
  #
  # Create the returned objects
  #

  # Estimated coefficient vector
  #coefficients <- unlist(lapply(results, function(x) x$coefficients))
  
  # Estimate variance-covariance matrix. This is block diagonal beacuse 
  # the covariances are available only for coefficients that are from the same
  # equation because 2SLS is not a system estimator
  
  #coefCov <- lavaan::lav_matrix_bdiag(lapply(results, function(x) x$coefCov))
  
  return(d)
  
  # The code below is the old 2SLS estimator code
  
  #-------------------------------------------------------#  
  # Loop over the equations and calculate 2sls estimates
  #-------------------------------------------------------#
  
  dvLabels    <- unlist(lapply(d,"[[","DVobs"))
  numEq       <- length(d)      # number of equations
  numCoef_i   <- numeric(numEq) # vector of number of coefficients in each eq.
  coefNames_i <- list() # any names of coefficients of each equation
  obsNames_i  <- list() # any names of observations of each equation
  terms_zy_i  <- list() # list of terms for Z and Y objects in each equation
  terms_v_i   <- list() # list of terms for V objects in each equation
  Y_i  <- list()        # list of dependent variable vectors in each eq.
  Z_i  <- list()        # list of explanatory variable matrices in each eq.
  V_i  <- list()        # list of instrumental variable matrices in each eq.
  browser()
  for(i in 1:numEq) { 
    
    data$formula <- NULL; evaldata <- NULL;
    data$formula   <- d[[i]]$EqFormula
    evaldata       <- eval(data)
    
    terms_zy_i[[i]] <- attr( evaldata, "terms" )
    Y_i[[i]]        <- stats::model.extract(evaldata, "response" )
    Z_i[[i]]        <- stats::model.matrix(terms_zy_i[[i]], evaldata)
    
    obsNames_i[[i]]  <- rownames(Z_i[[i]])
    
    numCoef_i[i]     <- ncol(Z_i[[i]])
    
    coefNames_i[[i]] <- paste0(d[[i]]$DVobs,"_", colnames(Z_i[[i]]))
    
    data$formula <- NULL; evaldata <- NULL;
    data$formula   <- d[[i]]$MIIVsFormula
    evaldata       <- eval(data)
    
    terms_v_i[[i]] <- attr(evaldata, "terms")
    V_i[[i]] <- model.matrix(terms_v_i[[i]], evaldata)
  }
  
  numCoef      <- sum(numCoef_i)
  coefNames    <- unlist(coefNames_i)
  
  #-------------------------------------------------------#  
  # Process any missing values.
  #-------------------------------------------------------#
  noNA        <- processMissing(Y_i, Z_i, V_i)
  numObsPerEq <- nrow(noNA)
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # For now remove any missing values.
  #-------------------------------------------------------#
  for(i in 1:numEq) {
    Y_i[[i]]   <- Y_i[[i]][noNA[,i]]
    attrAssign <- attributes(Z_i[[i]])$assign
    Z_i[[i]]   <- Z_i[[i]][noNA[, i], , drop = FALSE]
    attributes(Z_i[[i]])$assign <- attrAssign
    V_i[[i]] <- V_i[[i]][noNA[, i], , drop = FALSE ]
  }
  
  validObsNames_i <- lapply(Y_i,names)
  numObsEq_i      <- unlist(lapply(Y_i,length))
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Calculate stage 1 fitted regression matrices using 
  # raw data with the missing values removed.  
  #-------------------------------------------------------#
  Zhat_i <- calcStage1Fitted(Z_i, V_i)
  #-------------------------------------------------------#

  
  #-------------------------------------------------------#  
  # Calculate degrees of freedom.
  #-------------------------------------------------------#  
  df_list       <- calcDegreesOfFreedom(d, restrictions, numObsEq_i, numCoef_i)
  df_i          <- df_list$df
  numCoefRes    <- df_list$numCoefRes   # nCoefLiAll
  numCoefRes_i  <- df_list$numCoefRes_i # nCoefLiEq
  #-------------------------------------------------------# 
  
  #-------------------------------------------------------#  
  # Calculate stage 2 regression.
  #-------------------------------------------------------#
  b <- calcStage2coef(Y_i, Zhat_i, V_i, R, q, coefNames)
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Calculate the residual covariance matrix.
  #-------------------------------------------------------#  
  Omega <- calcResidCovMat(resids_i,numEq, noNA = noNA) 
  #-------------------------------------------------------# 
  
  #-------------------------------------------------------#  
  # Calculate coefficient covariance matrix.
  #-------------------------------------------------------#  
  coefCov <- calcCoefCovMat(Zhat_i, Omega, R = R, q = q, numEq, noNA, coefNames)
  #-------------------------------------------------------#  
  
}


#' Two-stage least square estimator for a system of equations
#' 
#' @param d a system of MIIV equation
#' @param sample.cov sample covariance matrix
#' @param sample.nobs number of observations
#' @param se should variance covariance matrix of the estimates be calculated
#' @param restrictions any equality constraints to be used in estimation

#'@keywords internal

miive.2sls.system <- function(d, sample.cov, sample.mean, sample.nobs, se, restrictions){
  
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
  
  names(coef) <- unlist(lapply(d, function(x) paste0(x$DVobs,"~", c("1", x$IVobs))))
  
  # Add coefficients to equations list.
  coefIndex <- unlist(lapply(seq_along(d), function(x) rep(x,(length(d[[x]]$IVobs)+1))))
  coefList  <- split(coef, coefIndex); names(coefList) <- rep("coefficients",length(d))
  d         <- lapply(seq_along(d), function(x) append(d[[x]], coefList[x])) 
  
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
    
    if (is.null(restrictions)){
      
      coefCov <- solve(XX1 %*% t(solve(sig)))
      
    } else {
      
      R0 <- matrix(0, ncol=nrow(R), nrow=nrow(R))
      coefCov <- solve(rbind(cbind((XX1 %*% t(solve(sig))), t(R)), 
                             cbind(R, R0)))[1:nrow(XX1), 1:nrow(XX1)]
    }
    
    d$coefCov <- coefCov
    
    
  }
  
}
