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
#' @param var.cov If \code{TRUE}, variance and covariance parameters are estimated.
#' @param se If "standard", conventional closed form standard errors are computed. If "boot" or "bootstrap", bootstrap standard errors are computed using standard bootstrapping.
#' @param bootstrap Number of bootstrap draws, if bootstrapping is used.
#' @param miivs.check Options to turn off check for user-supplied instruments validity as MIIVs.
#' @param ordered A vector of variable names to be treated as ordered factors in generating the polychoric correlation matrix.
#' @details 
#' 
#' \itemize{
#' \item{\code{instruments}} {
#'   Using the \code{instruments} option you can specify the MIIVs directly 
#'   for each equation in the model. To utilize this option you must first 
#'   define a list of instruments using the syntax displayed below. After 
#'   the list is defined, set the \code{instruments} argument equal to 
#'   the name of the list of MIIVs. Note, \code{instruments} are specified 
#'   for an equation, and not for a specific endogenous variable.
#'   }
#'\item{\code{se}} {
#'   When \code{se} is set to \code{"boot"} or \code{"bootstrap"} standard errors
#'   are computed using the pairs bootstrap. These standard errors are based
#'   on the standard deviation of successful bootstrap replications.  The 
#'   \code{z-value} and \code{P(>|z|)} assume the ratio of the coefficient 
#'   estimate to the bootstrap standard deviation approximates a normal 
#'   distribution.  Note, the Sargan test statistic is calculated from the 
#'   original sample and is not a bootstrap estimate.
#'   }
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
                  var.cov = FALSE, miiv.check = TRUE, ordered = NULL){
  
  #-------------------------------------------------------# 
  # A few basic sanity checks for user-supplied covariance 
  # matrices and mean vectors. Also, rescale cov.matrix if
  # cov.rescale set to TRUE.
  #-------------------------------------------------------# 
  if(is.null(data)){
    
    if (!is.null(ordered)){
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
  # remove equations that are not identified
  #-------------------------------------------------------#
  underid   <- sapply(d, function(eq) {
    length(eq$MIIVs) < length(eq$IVobs)
  })
  d.un <- d[ underid]; d   <- d[!underid]
  
  #-------------------------------------------------------# 
  # Remove variables from data that are not in model 
  # syntax and preserve the original column ordering.
  #-------------------------------------------------------# 
  if(!is.null(data)){
    data     <- data[,colnames(data) %in% unique(unlist(lapply(d, function(eq){
      c(eq$DVobs,eq$IVobs,eq$MIIVs)
    })))]
    data <- as.data.frame(data)
  }


  #-------------------------------------------------------# 
  # Process data. See documentation of processRawData. 
  #-------------------------------------------------------# 
  g <- processData(data, sample.cov, sample.mean, sample.nobs, ordered, pt)
  
  
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
    eq$exogenous   <- ifelse(any(g$var.exogenous[c(eq$DVobs, eq$IVobs, eq$MIIVs)]), TRUE, FALSE)
    eq$restricted  <- ifelse(r$eq.restricted[eq$DVobs], TRUE, FALSE)
    # For now throw an error if a categorical equatio contains an exogenous variable. 
    if (eq$categorical & eq$exogenous){
      stop(paste("miive: exogenous variables in",
                 "equations containing categorical variables (including MIIVs)",
                 "are not currently supported."))
    }
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
      "2SLS" = miive.2sls(d, d.un, g, r, est.only),
      "GMM"  = miive.gmm(d, d.un, g, r, est.only), # Not implemented
      # In other cases, raise an error
      stop(paste("Invalid estimator:", estimator, "Valid estimators are: 2SLS, GMM"))
  )
  
  #-------------------------------------------------------#
  # Estimate variance and covariance point
  #-------------------------------------------------------#
  if (var.cov){
    
    vcov.model <- createModelSyntax(results$eqn, pt)
    v <- estVarCovar(data, vcov.model, ordered)
    
  } else {
    
    v <- NULL
    
  }

  #-------------------------------------------------------#
  # Boostrap and substitute closed form SEs with boostrap SEs
  #-------------------------------------------------------#
  
  if(se == "boot" | se == "bootstrap"){
    
    # Do this outside of the bootstrap loop
    if (var.cov){
      
      if(length(d.un) > 0){
        stop(paste("MIIVsem: variance covariance esimtation not",
                   "allowed in the presence of underidentified",
                   "equations."))
      }
      
    }
    
    boot.results <- boot::boot(data, function(origData, indices){
      
      bsample <- origData[indices,]
      
      g <- processData(
        data = bsample, 
        sample.cov = NULL, 
        sample.mean = NULL, 
        sample.nobs = NULL, 
        ordered = ordered, 
        pt = pt
      )
      
      
      brep <- switch(
        estimator,
        "2SLS" = miive.2sls(d, d.un, g, r, est.only = TRUE),
        "GMM"  = miive.gmm(d, d.un, g, r, est.only = TRUE), # Not implemented
        # In other cases, raise an error
        stop(paste(
          "Invalid estimator:", 
          estimator, 
          "Valid estimators are: 2SLS, GMM")
        )
      )
      
      if (var.cov){
        
        num.vcov   <- nrow(pt[pt$op == "~~",])
        vcov.model <- createModelSyntax(brep$eqn, pt)
        vcov.coefs <- estVarCovar(bsample,vcov.model, ordered)$coefficients
        
        c(brep$coefficients, vcov.coefs)
        
      } else {
        
        brep$coefficients
        
      }
      
    }, bootstrap)
    
    # Replace the estimated variances of estimates 
    # with the boostrap estimates
    boot.mat       <- boot.results$t[complete.cases(boot.results$t),]
    results$bootstrap.true <- nrow(boot.mat)
    coefCov <- cov(boot.mat)
    rownames(coefCov) <- colnames(coefCov) <- c(
      names(results$coefficients), if (var.cov) names(v$coefficients) else NULL
    )
    results$coefCov <- coefCov

    # Store the boot object as a part of the result object. 
    # This is useful for calculating CIs or
    # other bootstrap postestimation.
    results$eqn <- lapply(results$eqn, function(eq) {
      if (eq$categorical) {
        eq$coefCov <- NA
      } else {
        eq$coefCov <- coefCov[
          paste0(eq$DVlat,"~",c("1", eq$IVlat)),
          paste0(eq$DVlat,"~",c("1", eq$IVlat)), 
          drop = FALSE
          ]
      }
      eq
    })

    results$boot <- boot.results
  }
  

  
  # assemble return object
  results$eqn.unid       <- d.un
  results$estimator      <- estimator
  results$se             <- se
  results$bootstrap      <- bootstrap
  results$call           <- match.call()
  results$r              <- r
  results$v              <- v
  
  class(results)  <- "miive"
  
  results
}

