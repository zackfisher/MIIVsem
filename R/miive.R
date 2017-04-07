#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using model-implied instrumental variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax or a \code{\link{miivs}} object. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by \code{as.data.frame} to data frame.
#' @param instruments A user-supplied list of valid MIIVs for each equation. See Example 2 below. 
#' @param sample.cov Numeric matrix. A sample variance-covariance matrix. The rownames and colnames must contain the observed variable names.
#' @param sample.mean A sample mean vector.
#' @param sample.nobs Number of observations in the full data frame.
#' @param sample.cov.rescale If \code{TRUE}, the sample covariance matrix provided by the user is internally rescaled by multiplying it with a factor (N-1)/N.
#' @param estimator Options \code{"2SLS"} or \code{"GMM"} for estimating the model parameters. Default is \code{"2SLS"}.
#' @param est.only If \code{TRUE}, only the coefficients are returned.
#' @param se If "standard", conventional closed form standard errors are computed. If "boot" or "bootstrap", bootstrap standard errors are computed using standard bootstrapping.
#' @param bootstrap Number of bootstrap draws, if bootstrapping is used.
#' @param miivs.check Options to turn off check for user-supplied instruments validity as MIIVs.
#' @param factor.vars A vector of variable names to be treated as ordered factors in generating the polychoric correlation matrix.
#' @details 
#' 
#' \itemize{
#' \item{\code{instruments}} {
#'   Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.
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
#' @seealso \link{MIIVsem}{miivs}
#'  
#' @export
miive <- function(model = model, data = NULL,  instruments = NULL,
                  sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL, 
                  sample.cov.rescale = TRUE, estimator = "2SLS", 
                  se = "standard", bootstrap = 1000L, est.only = FALSE, 
                  miiv.check=TRUE, factor.vars = NULL){
  
  #-------------------------------------------------------# 
  # A few basic sanity checks for user-supplied covariance 
  # matrices and mean vectors. Also, rescale cov.matrix if
  # cov.rescale set to TRUE.
  #-------------------------------------------------------# 
  if(is.null(data)){
    
    if (!is.null(factor.vars)){
      stop(paste("miive: raw data required when declaring factor variables."))
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
    
    if(sample.cov.rescale & !is.null(sample.cov)){
      sample.cov <- sample.cov * (sample.nobs-1)/sample.nobs
    }
  }
  
  #-------------------------------------------------------#  
  # Check class of model.
  #-------------------------------------------------------#
  if ( "miivs" == class(model) ){ 
    
    d  <- model$eqns 
    pt <- model$pt
    
  } else { 
    
    res <- miivs(model)
    d   <- res$eqns
    pt  <- res$pt
    
  } 
  
  #-------------------------------------------------------# 
  # parseInstrumentSyntax
  #-------------------------------------------------------#
  d  <- parseInstrumentSyntax(d, instruments, miiv.check)
  
  #-------------------------------------------------------# 
  # Remove variables from data that are not in model 
  # syntax and preserve the original column ordering.
  #-------------------------------------------------------# 
  if(!is.null(data)){
    data     <- data[,colnames(data) %in% unique(unlist(lapply(d, function(eq){
      c(eq$DVobs,eq$IVobs,eq$MIIVs)
    })))]
  }
  
  if(!is.null(factor.vars)){
    data[,factor.vars] <- lapply(data[,factor.vars], ordered)
  }
 

  
  #-------------------------------------------------------# 
  # Process data. See documentation of processRawData. 
  #-------------------------------------------------------# 
  g <- processData(data, sample.cov, sample.mean, sample.nobs, factor.vars, pt)
  
  
  #-------------------------------------------------------#  
  # Build Restriction Matrices.
  #-------------------------------------------------------#  
  r <- buildRestrictMat(d)
  
  #-------------------------------------------------------# 
  # Add some fields to d and check for any problematic
  # cases prior to estimation.  
  #-------------------------------------------------------# 
  d <- lapply(d, function (eq) {
    eq$missing <- ifelse(any(g$var.nobs[c(eq$DVobs, eq$IVobs, eq$MIIVs)] < g$sample.nobs), TRUE, FALSE)
    eq$categorical <- ifelse(any(g$var.categorical[c(eq$DVobs, eq$IVobs, eq$MIIVs)]), TRUE, FALSE)
    eq$restricted <- ifelse(r$eq.restricted[eq$DVobs], TRUE, FALSE)
    # throw an error if a categorical equation includes restrictions
    # if (eq$categorical & eq$restricted){
    #   stop(paste("miive: Restrictions on coefficients in",
    #              "equations containing categorical variables (including MIIVs)",
    #              "is not currently supported."))
    # }
    # if (eq$missing & eq$restricted){
    #   stop(paste("miive: Restrictions on coefficients in",
    #              "equations with missing values (including on the MIIVs)",
    #              "is not currently supported."))
    # }
    eq
  })
  
  #-------------------------------------------------------#
  # Generate results object.
  #-------------------------------------------------------#
  results <- switch(
    estimator,
      "2SLS" = miive.2sls(d, g, r, est.only),
      "GMM"  = miive.gmm(d, g, r, est.only), # Not implemented
      # In other cases, raise an error
      stop(paste("Invalid estimator:", estimator, "Valid estimators are: 2SLS, GMM"))
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
  results$estimator <- estimator
  class(results)  <- "miive"
  
  results
}

