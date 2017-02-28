#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using model-implied instrumental variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax or a \code{\link{miivs}} object. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by \code{as.data.frame} to data frame.
#' @param instruments A user-supplied list of valid MIIVs for each equation. See Example 2 below. 
#' @param sample.cov Numeric matrix. A sample variance-covariance matrix. The rownames and colnames must contain the observed variable names.
#' @param sample.mean A sample mean vector.
#' @param sample.nobs Number of observations if the full data frame is missing and only sample moments are given.
#' @param sample.cov.rescale If \code{TRUE}, the sample covariance matrix provided by the user is internally rescaled by multiplying it with a factor (N-1)/N.
#' @param estimator Options \code{"2SLS"}, \code{"GMM"} or \code{"PIV"} for estimating the model parameters. Default is \code{"2SLS"}.
#' @param control .
#' @param est.only If \code{TRUE}, only the coefficients are returned.
#' @param se If "standard", conventional closed form standard errors are computed. If "boot" or "bootstrap", bootstrap standard errors are computed using standard bootstrapping.
#' @param bootstrap Number of bootstrap draws, if bootstrapping is used.
#' @param piv.opts Options to pass to lavaan's \code{\link[lavCor]{lavCor}} function.
#' @param miivs.check Options to turn off check for user-supplied instruments validity as MIIVs.
#' @param factor.vars A vector of variable names to be treated as ordered factors in generating the polychoric correlation matrix.
#' @details 
#' 
#' \itemize{
#' \item{\code{instruments}} {Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.
#'   To utilize this option you must first define a list of instruments using the syntax displayed below. After the list is defined, 
#'   set the \code{instruments} argument equal to the name of the list of MIIVs. Note, \code{instruments} are specified for an equation,
#'   and not for a specific endogenous variable.}
#' }

#'
#' @return A list of class \code{miive} containing the following elements:
#'
#' \tabular{ll}{
#' \code{coefficients}\tab A named vector of parameter estimates\cr
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
#' @seealso \link{\code{miivs}}
#'  
#' @export
miive <- function(model = model, data = NULL,  instruments = NULL,
                  sample.cov = NULL, 
                  sample.mean = NULL, sample.nobs = NULL, 
                  sample.cov.rescale = TRUE, 
                  estimator = "2SLS", control = NULL, 
                  se = "standard", bootstrap = 1000L,
                  est.only = FALSE, piv.opts = c(estimator = "two.step", se = "standard"),
                  miivs.check=TRUE, factor.vars = NULL){
  
  #-------------------------------------------------------#  
  # Check class of model.
  #-------------------------------------------------------#
  if ( "miivs" == class(model) ){ 
    d  <- model$eqns 
    pt <- model$pt
  } 
  
  if ( "miivs" != class(model) ){ 
    res <- miivs(model)
    d   <- res$eqns
    pt  <- res$pt
  } 
  
  #-------------------------------------------------------# 
  # A few basic sanity checks for user-supplied covariance 
  # matrices and mean vectors. Also, rescale cov.matrix.
  #-------------------------------------------------------# 
  if(is.null(data)){
    
    if (!is.null(factor.vars)){
      stop(paste("miive: raw data required for using factor variables."))
    }
    
    if (!is.vector(sample.mean)){
      stop(paste("miive: sample.mean must be a vector."))
    }
    
    if (is.null(names(sample.mean))){
      stop(paste("miive: sample.mean vector must have names."))
    }
    
    if (!all.equal(names(sample.mean),colnames(sample.cov), check.attributes = FALSE)){
      stop(paste("miive: names of sample.mean vector and sample.cov matrix must match."))
    }
    
    # Rescale user-provided covariance matrix.
    if(sample.cov.rescale & ! is.null(sample.cov)){
      sample.cov <- sample.cov * (sample.nobs-1)/sample.nobs
    }
  }
  
  #-------------------------------------------------------# 
  # Check if any variables are categorical.
  #-------------------------------------------------------# 
  if (is.null(data)){
    factorIndex <- FALSE
  } else {
    factorIndex <- vapply(data,is.factor, c(is.factor=FALSE))
  }
  
  if(any(factorIndex)){
    
    # If factor.vars weren't given in the original input throw an error. 
    if (is.null(factor.vars)){
      stop(paste("miive: undeclared factors in dataset."))
    }
    
    # If factor.vars don't match those passed to miive by user throw an error. 
    if (!identical(sort(colnames(data)[factorIndex]),sort(factor.vars))){
      stop(paste("miive: undeclared factors found in dataset or factor.vars."))
    }

  }
  
  #-------------------------------------------------------# 
  # The estimation is done from covariance matrix and mean 
  # vector, so these are calculated first
  #-------------------------------------------------------# 
  if(!is.null(data)){
    data <- data[complete.cases(data),]
    sample.cov  <- cov(data)*(nrow(data)-1)/nrow(data)
    sample.nobs <- nrow(data)
    sample.mean <- colMeans(data)
  }
  
  #-------------------------------------------------------# 
  # parseInstrumentSyntax
  #-------------------------------------------------------# 
  #   input:  (1) miivs equation list 'd'
  #                returned from miivs search function.
  #           (2) instruments argument, null by default,
  #               if no instruments have been supplied.
  #           (3) miivs.check, default is true.  Can be 
  #               used to turn off miiv validity cehcks,
  #               for example, if instruments are not in
  #               model.
  #          
  #  return:  (1) updated 'd' 
  #
  #           The function contains checks to determine if (1)
  #           the user-supplied dependent variables exist in the
  #           set of estimating equations and (2) if the instruments
  #           provided by the user are valid MIIVs.
  #-------------------------------------------------------#
  d  <- parseInstrumentSyntax(d, instruments, miivs.check)
  
  
  #-------------------------------------------------------#  
  # Build Restriction Matrices.
  # returns NULL if there are no restrictions,
  # otherwise returns a list containing the 'R' matrix
  # and 'q' vector, as well as a vector 'cons' of 
  # the constrained coefficients.
  #-------------------------------------------------------#  
  restrictions <- buildRestrictMat(d)
  #-------------------------------------------------------#  
  
  # TODO: Check that all specified variables exists in the data and
  # throw an error if not.
  
  # TODO: Check if any variables are factors.  If so, 
  # set 'estimator' to "PIV". 
  
  #estimator <- ifelse(any(
  #   sapply(data[,unlist(unique(lapply(d,"[",c("DVobs","IVobs"))))], is.factor)
  #  ), "PIV", estimator)
    
  
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
                    "2SLS" = miive.2sls(d, data, sample.cov, sample.mean, sample.nobs, est.only, restrictions, piv.opts, factorIndex),
                    "GMM" = miive.gmm(d, data, sample.cov, sample.mean, sample.nobs, est.only, restrictions), # Not implemented
                    # In other cases, raise an error
                    stop(paste("Invalid estimator:", estimator,"Valid estimators are: 2SLS, GMM, PIV"))
  )
  
  #-------------------------------------------------------#
  # Estimate variance and covariance point
  #-------------------------------------------------------#
  
  # First fill the lavaan parTable with all regression style
  # coefficients and update the modely syntax.
  #modSyntax <- createModelSyntax(results$eqn, pt)

  # Obtain the variance and covariance point estimates.
  # if (estimator == "PIV"){
  #   varCoefs <- lavaan::parameterEstimates(lavaan::sem(
  #     modSyntax,
  #     data,
  #     ordered = colnames(data)[!apply(data,2,is.numeric)],
  #     estimator =  "ULS", # piv.opts["estimator"],
  #   ))
  # } else {
  #   varCoefs <- lavaan::parameterEstimates(lavaan::sem(
  #     modSyntax,
  #     sample.cov = sample.cov, 
  #     sample.nobs =sample.nobs,
  #     estimator =  "ULS"
  #   ))
  # }
  # varCoefs <- varCoefs[varCoefs$op == "~~" & !is.na(varCoefs$z),c("lhs", "rhs", "est")]
  # rownames(varCoefs) <- paste0(varCoefs$lhs, "~~",varCoefs$rhs)
  # results$varCoefs <- varCoefs
  results$varCoefs <- NULL
  #-------------------------------------------------------#
  # Boostrap and substitute closed form SEs with boostrap SEs
  #-------------------------------------------------------#
  
  if(se == "boot" | se == "bootstrap"){
    boot.results <- boot::boot(data,function(origData,indices){
      
      bsample <- origData[indices,]
      brep <- switch(estimator,
             "2SLS" = miive.2sls(d, bsample, sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL, est.only = TRUE, restrictions),
             "GMM" = miive.2sls(d, bsample, sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL, est.only = TRUE, restrictions)) # Not implemented
             # No need to raise an error here becaues we would have raised it earlier in any case
      brep$coefficients
    }, bootstrap)
    
    # Replace the estimated variances of estimates with the boostrap estimates
    results$coefCov <- cov(boot.results$t)
    
    # Store the boot object as a part of the result object. This is useful for calculating CIs or
    # other bootstrap postestimation.
    
    results$boot <- boot.results
  }
  
  # Keep the function call
  results$call <- match.call()
  
  class(results)  <- "miive"
  
  results
}

