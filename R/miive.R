#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using model-implied instrumental variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax or a \code{\link{miivs}} object. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by \code{as.data.frame} to data frame.
#' @param sample.cov Numeric matrix. A sample variance-covariance matrix. The rownames and colnames must contain the observed variable names.
#' @param sample.mean A sample mean vector.
#' @param sample.nobs Number of observations if the full data frame is missing and only sample moments are given.
#' @param sample.cov.rescale If \code{TRUE}, the sample covariance matrix provided by the user is internally rescaled by multiplying it with a factor (N-1)/N.
#' @param instruments A user-supplied list of valid MIIVs for each equation. See Example 2 below. 
#' @param estimator Options \code{"2SLS"} or \code{"GMM"} for estimating the model parameters. Default is \code{"2SLS"}.
#' @param control .
#' @param se If \code{TRUE}, standard errors are returned with estimates. 
#'
#' @return model
#' @return dat
#' @return modeeqns
#' 
#' @details 
#' \itemize{
#' \item{\code{instruments}} {Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.
#'   To utilize this option you must first define a list of instruments using the syntax displayed below. After the list is defined, 
#'   set the \code{instruments} argument equal to the name of the list of MIIVs. Note, \code{instruments} are specified for an equation,
#'   and not for a specific endogenous variable.}
#' }
#' 
#' @references 
#' 
#' Bollen, K. A. 1996.	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 
#' Bollen,	K. A. 2001.	Two-stage	Least	Squares	and	Latent	Variable	Models:	
#' Simultaneous	Estimation	and	Robustness	to	Misspecifications.
#' In	R.	Cudeck,	S.	Du	Toit,	and	D.	Sorbom	(Eds.),	Structural	
#' Equation	Modeling:	Present	and	Future,	A	Festschrift	in	Honor	of	Karl	
#' Joreskog	(pp. 119-138).	Lincoln,	IL: Scientific	Software.
#' 	
#'
#' @example example/bollen1989-miive1.R
#' @example example/bollen1989-miive2.R
#' @example example/bollen1989-miive3.R
#'  
#' @export
miive <- function(model = model, data = NULL, 
                  sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL, sample.cov.rescale = TRUE,
                  instruments = NULL, estimator = "2SLS", control = NULL, se = TRUE){
  
  #-------------------------------------------------------#  
  # Check class of model.
  #-------------------------------------------------------#
  if ( "miivs" == class(model) ){ d <- model} 
  if ( "miivs" != class(model) ){ d <- miivs(model)} 
  
  #-------------------------------------------------------# 
  # parseInstrumentSyntax
  #-------------------------------------------------------# 
  #   input:  (1) miivs equation list 'd'
  #                returned from miivs search function.
  #           (2) instruments argument, null by default,
  #               if no instruments have been supplied.
  #          
  #  return:  (1) updated 'd' 
  #
  #  details: 'd' updated with one new field:
  #           (1) MIIVsUsed is a vector of MIIVs to be used
  #           in estimation for equation i.
  #
  #           The function contains checks to determine if (1)
  #           the user-supplied dependent variables exist in the
  #           set of estimating equations and (2) if the instruments
  #           provided by the user are valid MIIVs.
  #-------------------------------------------------------#
  d  <- parseInstrumentSyntax(d, instruments)
  
  #-------------------------------------------------------#
  
  #-------------------------------------------------------# 
  # generateFormulas
  #-------------------------------------------------------# 
  #   input:  (1) miivs equation list 'd'
  #                returned from miivs search function.
  #          
  #  return:  (1) updated 'd' 
  #
  #  details: 'd' updated with two new fields:
  #           (1) EqFormula is two-sided formula for equation i.
  #           (2) MIIVsFormula is a one-sided formula
  #           characterizing the instruments for equation i.
  #
  #-------------------------------------------------------#
  d <- generateFormulas(d)
  
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Prepare data.  (REVISIT)
  #   
  #  To-do: Prepare for missing data.
  #         Prepare for covariance matrices
  #-------------------------------------------------------#
  
  # Mikko: I do not think that this is really needed
  # different estimators may treat missing data differently,
  # so I would defer missing data treatment to actual estimation
  # see my comments in the 2SLS estimator
  
  # mf <- prepareRawData(data)
  
  if(sample.cov.rescale & ! is.null(sample.cov)){
    sample.cov <- sample.cov * (sample.nobs-1)/sample.nobs
  }
  
  
  #-------------------------------------------------------#  
  # Build Restriction Matrices.
  # returns NULL if there are no restrictions,
  # otherwise returns a list containing the 'R' matrix
  # and 'q' vector, as well as a vector 'cons' of 
  # the constrained coefficients.
  #-------------------------------------------------------#  
  restrictions <- buildRestrictMat(d)
  
  if (is.null(restrictions)){ 
    
    R <- NULL; q <- NULL
    
  } else { 
    
    R <- restrictions$R; q <- restrictions$q 
    
  }
  #-------------------------------------------------------#  
  
  
  
  #-------------------------------------------------------#
  #
  # MIIV estimation using estimation functions. An 
  # estimation function retuns a list containing the
  # following elements:
  #
  # coefficients - a vector of estimated coefficients
  # coefCov - variance covariance matrix of the estimates
  #           (optional, depending on the se argument)
  # eqn - a list of miiv equations and estimation results, including
  #   Z - a matrix projecting instruments to the instrumented
  #     variables.
  #   
  #
  # possibly other elements
  # 
  # 
  #-------------------------------------------------------#
  
  results <- switch(estimator,
                    "2SLS" = miive.2sls(d, data, sample.cov, sample.mean, sample.nobs, se, restrictions),
                    "GMM" = miive.gmm(d, data, sample.cov, sample.mean, sample.nobs, se, restrictions), # Not implemented
                    # In other cases, raise an error
                    stop(paste("Invalid estimator:", estimator,"Valid estimators are: 2SLS, GMM"))
                    )

  #-------------------------------------------------------#  
  #
  # Calculate additional statistics that are common for
  # all estimators. These are calculated for each equation
  #-------------------------------------------------------#
  
  # Observation level statistics are calculated only if
  # raw data are available
  
  if(!is.null(data)){
    
    results$fitted <- do.call(cbind, lapply(results$eqn, function(eq){
      
      # Zack: Fitted values are obtained by multiplying the observed endogenous 
      # variables by the second stage coefficients.
      fitted <- drop(cbind("1"=1, data[,eq$IVobs]) %*% eq$coefficients)
      
      colnames(fitted) <-eq$DVlat
      
      return(fitted)
      
      # Fitted values are obtained by multiplying the observed
      # variables with the first stage coefficients and then
      # the second stage coefficients
      
      #fitted <- cbind("1"=1, designMatrix[,colnames(eq$Z)[-1]] %*% 
      #          t(eq$Z[-1,-1, drop = FALSE])) %*% eq$coefficients
      
      #colnames(fitted) <-eq$DVlat
      
      #fitted
    }))
    
  results$resCov <-   Sigma <- lavaan::lav_matrix_bdiag(lapply(d, function(eq){
    data[,sapply(results$eqn,function(eq){eq$DVobs})]   })) - results$fitted
  
  colnames(results$residuals) <- colnames(results$fitted)

  }
  
  # Keep the function call
  results$call <- match.call()
  
  class(results)  <- "miive"
  
  results
}

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