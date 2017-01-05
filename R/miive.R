#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using model-implied instrumental variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax or a \code{\link{miivs}} object. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by \code{as.data.frame} to data frame.
#' @param sample.cov Numeric matrix. A sample variance-covariance matrix. The rownames and colnames must contain the observed variable names.
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
miive <- function(model = model, data = NULL, sample.cov = NULL, instruments = NULL, 
                  estimator = "2SLS", control = NULL, se = TRUE){
  
  #-------------------------------------------------------#  
  # Check class of model.
  #-------------------------------------------------------#
  
  # Mikko: Would it not be simpler to just refer to model as model through the
  # code?
  
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
  
  d     <- parseInstrumentSyntax(d, instruments)
  
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
  
  d     <- generateFormulas(d)
  
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Prepare data.  (REVISIT)
  #   
  #  To-do: Prepare for missing data.
  #         Prepare for covariance matrices
  #-------------------------------------------------------#
  
  mf <- prepareRawData(data)
  
  #-------------------------------------------------------#
  #
  # MIIV estimation using estimation functions. An 
  # estimation function retuns a list containing the
  # following elements:
  #
  # coefficients - a vector of estimated coefficients
  # vcov - variance covariance matrix of the estimates
  #        (optional, depending on the se argument)
  # possibly other elements
  # 
  # 
  #-------------------------------------------------------#
  
  results <- switch(estimator,
                    "2SLS" = miive.2sls(d, mf, sample.cov, se),
                    "GMM" = miive.gmm(d, mf, sample.cov, se), # Not implemented
                    # In other cases, raise an error
                    stop(paste("Invalid estimator:", estimator,"Valid estimators are: 2SLS, GMM"))
                    )

  #-------------------------------------------------------#  
  #
  # Calculate additional statistics that are common for
  # all estimators
  #
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Calculate the residuals.
  #-------------------------------------------------------# 
  resids_i <- calcResiduals(b, Y_i, Z_i, numCoef_i)
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
  
  #-------------------------------------------------------#  
  # Fitted.values
  #-------------------------------------------------------#  
  fitted.values <- bdiag(Z_i) %*% b  
  #-------------------------------------------------------#  
  
  #-------------------------------------------------------#  
  # Organize results and calculate additional results
  #-------------------------------------------------------#
  b_i <- split(b, rep(1:length(numCoef_i), numCoef_i))
  results <- list()
  
  for(i in 1:numEq) {
    results$eq[[i]]          <- list()
    results$eq[[i]]$eqnNum   <- i         # equation number
    results$eq[[i]]$eqnLabel <- NA        # eqnLabels[[i]]
    results$eq[[i]]$method   <- estimator
    
    # Residuals
    results$eq[[i]]$residuals        <- resids_i[[i]]
    names(results$eq[[i]]$residuals) <- validObsNames_i[[i]]
    
    # Coefficients
    results$eq[[i]]$coefficients         <- b_i[[i]]
    #names(results$eq[[i]]$coefficients)  <- coefNames_i[[i]]
    
    # Coefficient Covariance Matrices
    results$eq[[i]]$coefCov  <- as.matrix(
      coefCov[(1+sum(numCoef_i[1:i])-numCoef_i[i]):(sum(numCoef_i[1:i])),
              (1+sum(numCoef_i[1:i])-numCoef_i[i]):(sum(numCoef_i[1:i]))] 
    )
    colnames(results$eq[[i]]$coefCov)    <- coefNames_i[[i]]
    rownames(results$eq[[i]]$coefCov)    <- coefNames_i[[i]]
    
    # Fitted Values
    results$eq[[i]]$fitted.values <- rep(NA, numObsPerEq)
    results$eq[[i]]$fitted.values[noNA[,i]] <-
      fitted.values[(1+sum(numObsEq_i[1:i])-numObsEq_i[i]):(sum(numObsEq_i[1:i]))]
    
    names(results$eq[[i]]$fitted.values) <- obsNames_i[[i]]
    
    # Terms
    results$eq[[i]]$terms <- terms_zy_i[[i]]
    
    # Rank
    results$eq[[i]]$rank  <- numCoefRes_i[i]
    
    # Total # of coef. in system
    results$eq[[i]]$nCoef.sys <- numCoef
    
    # Total # of rest. coef. in system
    results$eq[[i]]$rank.sys  <-  numCoefRes    # nCoefLiAll
    
    # rank = number of linear independent coefficients of the entire system
    results$eq[[i]]$df.residual  <- df_i[i]    
    results$eq[[ i ]]$df.residual.sys  <- sum(numObsEq_i) - numCoefRes
    results$eq[[i]]$inst      <- d[[i]]$MIIVsUsed
    results$eq[[i]]$termsInst <- terms_v_i[[i]]
    class(results$eq[[i]]) <- "miive.equation"
  }

  # System coefficients
  results$coefficients <- as.numeric(drop(b))
  names(results$coefficients) <- coefNames
  
  # Full Coefficient Covariance Matrix
  results$coefCov <- as.matrix(coefCov)
  dimnames(results$coefCov) <- list(coefNames, coefNames)
  
  # Residual Covariance Matrix 
  results$residCovEst <- as.matrix(Omega)
  dimnames(results$residCovEst) <- list(dvLabels, dvLabels)

  # Residual Covarance Matrix
  results$residCov <- function(resids_i, numEq, noNA, oneSigma = FALSE,
                               centered = FALSE, diagOnly = FALSE)
  dimnames(results$residCov) <- list(dvLabels, dvLabels)
  
  results$method  <- estimator
  results$rank    <- numCoefRes
  # rank = total number of linear independent coefficients of all equations
  results$df.residual <- sum(numObsEq_i) - numCoefRes
  # degrees of freedom of the whole system
  results$restrict.matrix <- R
  results$restrict.rhs    <- q
  results$control <- control
  class(results)  <- "miive"
  
  results
}
