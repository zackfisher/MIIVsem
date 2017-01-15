#'@method summary miive
#'@export
summary.miive <- function(x,..){
  
  # Mikko: I do not think that any of this is really needed in miive
  # because all the information is available in the coefficient vector,
  # model specification and the variance-covariance matrix
  # of the estimates. If we want estimation results as a list
  # of equations, that could be a part of a summary method
  # Note that this code does not currently work.  
  
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
  # Organize results and calculate additional results
  #-------------------------------------------------------#
  b_i <- split(b, rep(1:length(numCoef_i), numCoef_i))
  
  for(i in 1:numEq) {
    x$eq[[i]]          <- list()
    x$eq[[i]]$eqnNum   <- i         # equation number
    x$eq[[i]]$eqnLabel <- NA        # eqnLabels[[i]]
    x$eq[[i]]$method   <- estimator
    
    # Residuals
    x$eq[[i]]$residuals        <- resids_i[[i]]
    names(x$eq[[i]]$residuals) <- validObsNames_i[[i]]
    
    # Coefficients
    x$eq[[i]]$coefficients         <- b_i[[i]]
    #names(x$eq[[i]]$coefficients)  <- coefNames_i[[i]]
    
    # Coefficient Covariance Matrices
    x$eq[[i]]$coefCov  <- as.matrix(
      coefCov[(1+sum(numCoef_i[1:i])-numCoef_i[i]):(sum(numCoef_i[1:i])),
              (1+sum(numCoef_i[1:i])-numCoef_i[i]):(sum(numCoef_i[1:i]))] 
    )
    colnames(x$eq[[i]]$coefCov)    <- coefNames_i[[i]]
    rownames(x$eq[[i]]$coefCov)    <- coefNames_i[[i]]
    
    # Fitted Values
    x$eq[[i]]$fitted.values <- rep(NA, numObsPerEq)
    x$eq[[i]]$fitted.values[noNA[,i]] <-
      fitted.values[(1+sum(numObsEq_i[1:i])-numObsEq_i[i]):(sum(numObsEq_i[1:i]))]
    
    names(x$eq[[i]]$fitted.values) <- obsNames_i[[i]]
    
    # Terms
    x$eq[[i]]$terms <- terms_zy_i[[i]]
    
    # Rank
    x$eq[[i]]$rank  <- numCoefRes_i[i]
    
    # Total # of coef. in system
    x$eq[[i]]$nCoef.sys <- numCoef
    
    # Total # of rest. coef. in system
    x$eq[[i]]$rank.sys  <-  numCoefRes    # nCoefLiAll
    
    # rank = number of linear independent coefficients of the entire system
    x$eq[[i]]$df.residual  <- df_i[i]    
    x$eq[[ i ]]$df.residual.sys  <- sum(numObsEq_i) - numCoefRes
    x$eq[[i]]$inst      <- d[[i]]$MIIVsUsed
    x$eq[[i]]$termsInst <- terms_v_i[[i]]
    class(x$eq[[i]]) <- "miive.equation"
  }
  
  # System coefficients
  x$coefficients <- as.numeric(drop(b))
  names(x$coefficients) <- coefNames
  
  # Full Coefficient Covariance Matrix
  x$coefCov <- as.matrix(coefCov)
  dimnames(x$coefCov) <- list(coefNames, coefNames)
  
  # Residual Covariance Matrix 
  x$residCovEst <- as.matrix(Omega)
  dimnames(x$residCovEst) <- list(dvLabels, dvLabels)
  
  # Residual Covarance Matrix
  x$residCov <- function(resids_i, numEq, noNA, oneSigma = FALSE,
                         centered = FALSE, diagOnly = FALSE)
    dimnames(x$residCov) <- list(dvLabels, dvLabels)
  
  x$method  <- estimator
  x$rank    <- numCoefRes
  # rank = total number of linear independent coefficients of all equations
  x$df.residual <- sum(numObsEq_i) - numCoefRes
  # degrees of freedom of the whole system
  x$restrict.matrix <- R
  x$restrict.rhs    <- q
  x$control <- control
  
  return(x)
}