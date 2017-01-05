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

miive.2sls <- function(d, data, sample.cov, sample.mean, sample.nobs, se){
  
  # Estimate each equation
  
  results <- lapply(d, function(eq){
    
    # If we have data, then calculate a sample covariance matrix for this equation
    
    if(!is.null(data)){
      # Which variables are used
      v <- c(eq$DVobs, eq$IVobs, eq$MIIVsUsed)
      
      # Select the variables and use only complete cases
      s <- data[,v]
      s <- s[complete.cases(s),]
      
      miive.2sls.oneEq(eq, cov(s)*(nrow(s)-1)/nrow(s), colMeans(s), nrow(s), se)
    }
    
    # Else use the sample covariance matrix given as argument
    
    else{
      miive.2sls.oneEq(eq, sample.cov, sample.mean, sample.nobs, se)
    }    
  })  
  
  #
  # Create the returned objects
  #

  # Estimated coefficient vector
  coefficients <- unlist(lapply(results, function(x) x$coefficients))
  
  # Estimate variance-covariance matrix. This is block diagonal beacuse 
  # the covariances are available only for coefficients that are from the same
  # equation because 2SLS is not a system estimator
  
  coefCov <- lavaan::lav_matrix_bdiag(lapply(results, function(x) x$coefCov))
  
  return(list(coefficients = coefficients,
              coefCov = coefCov, 
              eqn = results # return the equation level results for convenience
              )
         )
  
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
  # Build Restriction Matrices.
  # returns NULL if there are no restrictions,
  # otherwise returns a list containing the 'R' matrix
  # and 'q' vector, as well as a vector 'cons' of 
  # the constrained coefficients.
  #-------------------------------------------------------#  
  restrictions <- buildRestrictMat(d)
  if (is.null(restrictions)){ R <- NULL; q <- NULL
  } else { R <- restrictions$R; q <- restrictions$q }
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


#' Two-stage least square estimator for a single equation
#' 
#' @param eq a MIIV equation
#' @param sample.cov sample covariance matrix
#' @param se should variance covariance matrix of the estimates be calculated
#' @param sample.nobs number of observations

#'@keywords internal

miive.2sls.oneEq <- function(eq, sample.cov, sample.mean, sample.nobs, se){
  
  # Stage 1
  
  # Crossproducts needed for OLS. These are not the actual crossproduct,
  # but are are divided by sample.nobs this does not influence the estimates
  # because all crossproducts are scaled similarly
  
  XX1 <- sample.cov[eq$MIIVs,eq$MIIVs] + sample.mean[eq$MIIVs] %o% sample.mean[eq$MIIVs]
  XX1 <- rbind(c("1"=1,sample.mean[eq$MIIVs]),cbind("1"=sample.mean[eq$MIIVs],XX1))

  XY1 <- rbind(sample.mean[eq$IVobs],sample.cov[eq$MIIVs,eq$IVobs] + sample.mean[eq$MIIVs]%o%sample.mean[eq$IVobs])
  
  coef1 <-t(solve(XX1, XY1))
  
  # Create a projection matrix Z projeting the MIIVs as fitted values
  # The latent DV is just replaced with the marker indicator
  
  Z <- rbind(0,cbind(0,coef1))
  Z[1,1] <- 1          
  
  colnames(Z)[1] <- eq$DVobs
  rownames(Z) <- c(eq$DVlat, eq$IVlat)

  # Stage 2
  
  # Calculate new crossproduct matrices for the second stage

  XX2 <- Z[-1,-1, drop = FALSE] %*% XX1 %*% t(Z[-1,-1, drop = FALSE])
  cp <- as.vector(coef1 %*% c(1,sample.mean[eq$MIIVs]))
  XX2 <- rbind(c("1"=1,cp),
               cbind("1"=cp,XX2))
  XY2 <- c(sample.mean[eq$DVobs],rowSums(Z[-1,-1, drop = FALSE] %*% XY1))
  
  coef2 <- solve(XX2,XY2)
  names(coef2) <- paste(eq$DVlat, c("1",eq$IVlat),sep= "~")
  
  # Add the results to the end of the MIIV equations
  eq$coefficients <- coef2
  eq$Z <- Z
  
  # Variance covariance matrix if requested
  
  if(se){
    
    # Error variance estimate. We eliminate the intercept from the coefficients
    # and crossproducts when calculating this

    sigma2 <- as.vector(sample.cov[eq$DVobs,eq$DVobs] - coef2[-1] %*% XX2[-1,-1] %*% coef2[-1])

    eq$coefCov <- sigma2 * solve(XX2*sample.nobs)
    
    colnames(eq$coefCov) <- rownames(eq$coefCov) <- names(coef2)
  }
  
  eq
}
