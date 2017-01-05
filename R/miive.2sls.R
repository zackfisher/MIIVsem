#' Two-stage least square estimator
#'
#'@keywords internal

miive.2sls <- function(d, mf, sample.cov, se){
  
  #-------------------------------------------------------#  
  # Prepare lists for terms and matrices.
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
  
  for(i in 1:numEq) { 
    
    mf$formula <- NULL; evalMF <- NULL;
    mf$formula   <- d[[i]]$EqFormula
    evalMF       <- eval(mf)
    
    terms_zy_i[[i]] <- attr( evalMF, "terms" )
    Y_i[[i]]        <- stats::model.extract(evalMF, "response" )
    Z_i[[i]]        <- stats::model.matrix(terms_zy_i[[i]], evalMF)
    
    obsNames_i[[i]]  <- rownames(Z_i[[i]])
    
    numCoef_i[i]     <- ncol(Z_i[[i]])
    
    coefNames_i[[i]] <- paste0(d[[i]]$DVobs,"_", colnames(Z_i[[i]]))
    
    mf$formula <- NULL; evalMF <- NULL;
    mf$formula   <- d[[i]]$MIIVsFormula
    evalMF       <- eval(mf)
    
    terms_v_i[[i]] <- attr(evalMF, "terms")
    V_i[[i]] <- model.matrix(terms_v_i[[i]], evalMF)
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
  b <- calcStage2Coefs(Y_i, Zhat_i, V_i, R, q, coefNames)
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
  
  return(list(coefficients = b))
}