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
#' @param est.only If \code{TRUE}, only the coefficients are returned.
#'
#' @details 
#' 
#' \itemize{
#' \item{\code{coefficients}} {Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.
#'   To utilize this option you must first define a list of instruments using the syntax displayed below. After the list is defined, 
#'   set the \code{instruments} argument equal to the name of the list of MIIVs. Note, \code{instruments} are specified for an equation,
#'   and not for a specific endogenous variable.}
#' }

#'
#' @return A list of class \code{miive} containing the following elements:
#'
#' \tabular{ll}{
#' \code{coefficients}\tab A named vector of parameter estimats\cr
#' \code{coefCov}\tab A variance-covariance matrix of the parameter estimates\cr
#' \code{residCov}\tab A residual variance-covariance matrix\cr
#' \code{eqn}\tab Equation level estimation resutls and statistics\cr
#' \code{call}\tab The matched call\cr
#'}
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
#' @seealso \link{code{miivs}}
#'  
#' @export
miive <- function(model = model, data = NULL, sample.cov = NULL, 
                  sample.mean = NULL, sample.nobs = NULL, 
                  sample.cov.rescale = TRUE, instruments = NULL, 
                  estimator = "2SLS", control = NULL, est.only = FALSE){
  
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
  #           The function contains checks to determine if (1)
  #           the user-supplied dependent variables exist in the
  #           set of estimating equations and (2) if the instruments
  #           provided by the user are valid MIIVs.
  #-------------------------------------------------------#
  d  <- parseInstrumentSyntax(d, instruments)
  
  
  #-------------------------------------------------------#  
  # Rescale user-provided covariance matrix.
  #-------------------------------------------------------#
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
  #-------------------------------------------------------#  
  
  #-------------------------------------------------------#
  # estimator: miive.2sls()
  #-------------------------------------------------------#
  # MIIV estimation using estimation functions. An 
  # estimation function retuns a list containing the
  # following elements:
  #
  # coefficients - a vector of estimated coefficients
  # coefCov      - variance covariance matrix of the estimates
  #                (optional, depending on the est.only argument)
  # residCov     - variance covariance matrix of the equation
  #                disturbances (optional, depending on the 
  #                est.only argument).
  # eqn -          a list of miiv equations and equation level
  #                estimation results. fields added include:
  #                  coefficients: coefficients
  #                  sigma       : equation sigma^2
  #                  sargan      : Sargan's test statistic
  #                  sargan.df   : df freedom for Sargan's
  #-------------------------------------------------------#
  
  results <- switch(estimator,
                    "2SLS" = miive.2sls(d, data, sample.cov, sample.mean, sample.nobs, est.only, restrictions),
                    "GMM" = miive.gmm(d, data, sample.cov, sample.mean, sample.nobs, est.only, restrictions), # Not implemented
                    # In other cases, raise an error
                    stop(paste("Invalid estimator:", estimator,"Valid estimators are: 2SLS, GMM"))
                    )

  # Keep the function call
  results$call <- match.call()
  
  class(results)  <- "miive"
  
  results
}

