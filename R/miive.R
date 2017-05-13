#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using model-implied instrumental variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax or a
#'        \code{\link{miivs}} object returned by the \code{\link{miivs}}
#'        function. See details for more information about permissible
#'        operators and example syntax.
#' @param data A data frame, list or environment or an object coercible 
#'        by \code{as.data.frame} to data frame.
#' @param instruments A user-supplied list of instruments for each 
#'        equation. See below and the \code{miivs.out} argument of 
#'        \code{\link{miivs}} for more information on the correct 
#'        input format. External (auxiliary) instruments can be 
#'        supplied, however, the \code{miiv.check} argument must
#'        be set to \code{FALSE}. 
#' @param sample.cov Numeric matrix. A sample variance-covariance 
#'        matrix. The rownames and colnames must contain each of
#'        the observed variable names indicated in the model syntax.
#' @param sample.mean A sample mean vector.
#' @param sample.nobs Number of observations in the full data frame.
#' @param sample.cov.rescale If \code{TRUE}, the sample covariance matrix
#'        provided by the user is internally rescaled by multiplying 
#'        it with a factor (N-1)/N.
#' @param estimator Options \code{"2SLS"} or \code{"GMM"} for estimating the
#'        model parameters. Default is \code{"2SLS"}. Currently, only 
#'        \code{2SLS} is supported.
#' @param est.only If \code{TRUE}, only the coefficients are returned.
#' @param var.cov If \code{TRUE}, variance and covariance parameters are
#'        estimated.
#' @param var.cov.estimator The estimator to use for variance and covariance
#'        parameters.
#' @param se If "standard", conventional closed form standard errors are
#'        computed. If "boot" or "bootstrap", bootstrap standard errors 
#'        are computed using standard bootstrapping.
#' @param bootstrap Number of bootstrap draws, if bootstrapping is used.
#' @param missing If "listwise"...
#' @param miiv.check Options to turn off check for user-supplied 
#'        instruments validity as MIIVs.
#' @param ordered A vector of variable names to be treated as ordered factors
#'        in generating the polychoric correlation matrix.
#' @param wald Logical indicating whether or not to perform a wald test on any
#'        linear restrictions.
#' 
#' @details 
#' 
#' \itemize{
#' \item{\code{model}} {
#' 
#'   A model specified using a subset of the uses a subset of the model syntax 
#'   employed by \pkg{lavaan} or a \code{miivs} object return by the 
#'   \code{miivs} functions. See the \code{model} argumen within the 
#'   \code{\link[lavaan]{lavaanify}} function for more information. 
#'   The following model syntax operators are currently supported: \code{=~},
#'   \code{~}, \code{~~} and \code{*}. See below for details 
#'   on default behavior descriptions of how to specify the scaling 
#'   indicator in latent variable models and impose equality constraints 
#'   on the parameter estimates. 
#'   
#'   \strong{Example using Syntax Operators}
#'   
#'   In the model below, 'L1 \code{=~} Z1 + Z2 + Z3'  indicates the 
#'   latent variable L1 is measured by 3 indicators, Z1, Z2, and Z3. Likewise,
#'   L2 is measured by 3 indicators, Z4, Z5, and Z6. The statement
#'   'L1 \code{~} L2' specifies latent  variable L1 is regressed on latent 
#'   variable L2. 'Z1 \code{~~} Z2' specifies the error of indicator 
#'   Z2 is allowed to covary with the error of indicator Z3. The label
#'   L3 appended to Z3 and Z6 in the measurement model equations 
#'   constrains the factor loadings for Z3 and Z6 to equality. Additional 
#'   details on constraints see Equality Constraints  and Parameter 
#'   Restrictions.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + Z2 + L3*Z3
#'      L2 =~ Z4 + Z5 + L3*Z6
#'      L1  ~ L2
#'      Z2 ~~ Z3
#'   '}
#'  
#'   \strong{Scaling Indicators}
#'   
#'   Following the \pkg{lavaan} model syntax, latent variables are defined 
#'   using the \code{=~} operator.  For first order factors, the scaling 
#'   indicator chosen is the first observed variable on the RHS of an 
#'   equation. For the model  below \code{Z1} would be chosen as the 
#'   scaling indicator for \code{L1} and \code{Z4} would be chosen as 
#'   the scaling indicator for \code{L2}.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + Z2 + Z3
#'      L2 =~ Z4 + Z5 + Z6
#'   '}
#'   
#'   \strong{Higher-order Factor Models}
#'   
#'   For example, in the model below, the  scaling indicator for the 
#'   higher-order factor \code{H1} is taken to be \code{Z1}, the scaling 
#'   indicator that would have been assigned to the first lower-order 
#'   factor \code{L1}.
#'   
#'   \preformatted{model <- '
#'      H1 =~ L1 + L2 
#'      L1 =~ Z1 + Z2 + Z3
#'      L2 =~ Z4 + Z5 + Z6
#'   '}
#'   
#'   \strong{Equality Constraints and Parameter Restrictions}
#'   
#'   Within- and across-equation equality constraints on the factor loading
#'   and regression coefficients can be imposed directly in the model syntax. 
#'   To specify equality constraints between different parameters equivalent
#'   labels should be prepended to the variable name using the 
#'   \code{*} operator. For example, we could constrain the factor 
#'   loadings for two non-scaling indicators of latent factor \code{L1} to 
#'   equality using the following  model syntax.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + B1*Z2 + B1*Z3
#'      L2 =~ Z4 + Z5 + Z6
#'   '}
#'   
#'   The factor loading and regression coefficients can also be constrained
#'   to specific numeric values in a similar fashion. Below we constrain
#'   the regression coefficient  of \code{L1} on \code{L2} to \code{1}.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + Z2 + Z3
#'      L2 =~ Z4 + Z5 + Z6
#'      L3 =~ Z7 + Z8 + Z9
#'      L1  ~ 1*L2 + L3
#'   '}
#'   
#'   \strong{Model Defaults}
#'   
#'   In addition to those relationships specified in the model syntax 
#'   \pkg{MIIVsem} will automatically include the intercepts of any 
#'   observed or latent variables. Covariances among exogenous latent
#'   and observed  variables are included by default. 
#'   Where appropriate the covariances of latent and observed dependent 
#'   variables are also automatically included in the model specification.
#'   These defaults correspond to those used by \pkg{lavaan} and 
#'   \code{auto = TRUE}. 
#'   
#'   For example, in the model below the errors of Z1 and Z4 are allowed to 
#'   covary, in addition to all pairs of exogenous variables.
#'   
#'   \preformatted{model <- '
#'     Z1 ~ Z2 + Z3 
#'     Z4 ~ Z5 + Z6
#'   '}
#'   
#'   \strong{Invalid Specifications}
#'   
#'   Certain model specifications are not currently supported.  For example,
#'   the scaling indicator of a latent variable is not permitted to
#'   cross-load on another latent variable. For example, in the model below
#'   \code{Z1}, the scaling indicator for L1, cross-loads on the latent 
#'   variable \code{L2}. Executing a search on the model below will 
#'   result in the warning: \emph{miivs: scaling indicators with a factor 
#'   complexity greater than 1 are not currently supported}.
#'   
#'   \preformatted{model <- '
#'     L1 =~ Z1 + Z2 + Z3
#'     L2 =~ Z4 + Z5 + Z6 + Z1
#'   '}
#'   
#'   }
#'   
#'   \item{\code{data}} {
#'   
#'   A data.frame containing all of the observed variables specified in the
#'   model syntax. 
#'   }
#'   
#'   \item{\code{instruments}} {
#' 
#'   Using the \code{instruments} option you can specify the MIIVs directly 
#'   for each equation in the model. To utilize this option you must first 
#'   define a list of instruments using the syntax displayed below. Here,
#'   the dependent variable for each equation is listed on the LHS of the
#'   \code{"~"} operator. In the case of latent variable equations, the
#'   dependent variable is the scaling indicator associated with that
#'   variable. The instruments are then given on the RHS, seperated
#'   by \code{"+"} signs. For example, 
#'   
#'   \preformatted{myInstruments <- '
#'      y1 ~ z1 + z2 + z3
#'      y2 ~ z4 + z5
#'   '
#'   }
#'     
#'   After this list is defined, set the \code{instruments} argument equal to 
#'   the name of the list of instruments (e.g. \code{myInstruments}). 
#'   Note, that \code{instruments} are specified for an equation, 
#'   and not for a specific endogenous variable. If only a subset of dependent
#'   variables are listed in the instruments argument, only those  equations 
#'   will be estimated.  If external or auxilliary instruments (instruments 
#'   not otherwise included in the model) the \code{miiv.check} argument 
#'   should be set to \code{FALSE}.
#'   }
#'   
#'   \item{\code{sample.cov}} {
#'   
#'   The user may provide a sample covariance matrix in lieu of raw data.
#'   If \code{sample.cov} is not \code{NULL} the user must also supply a
#'   vector of sample means (\code{sample.mean}), and the number of sample
#'   observations (\code{sample.nobs}) from which the means and covariances 
#'   were calculated.  Currently, \pkg{MIIVsem} does not support bootstrap 
#'   standard errors or polychoric instrumental variable esimtation when
#'   the sample moments, rather than raw data, are used as input.  
#'   }
#'   
#'   \item{\code{sample.mean}} {
#'   
#'   A named vector of length corresponding to the row and column dimensions
#'   of the \code{sample.cov} matrix.  Names must also match those in the
#'   \code{sample.cov}.  
#'   }
#'   
#'   \item{\code{sample.cov.rescale}} {
#'   
#'   Default is \code{TRUE}.  If the sample covariance matrix provided 
#'   by the user should be internally rescaled by multiplying it with a 
#'   factor (N-1)/N.
#'   }
#'   
#'   \item{\code{estimator}} {
#'   
#'   The default estimator is \code{2SLS}. For equations with continuous
#'   variables only and no restrictions the estimatates are identical to 
#'   those described in Bollen (1996, 2001). If restrictions are present 
#'   a restricted MIIV-2SLS estimator is implemented using methods 
#'   similar to those described by Greene (2000) but adapted for 
#'   covariance based estimation. 2SLS coefficients and overidentifcation
#'   tests are constructed using means and covariances only for 
#'   increased computational efficiency. 
#'   
#'     If an equation contains ordered 
#'   categorical variables, declared in the \code{ordered} argument,
#'   the PIV estimator described by Bollen and Maydeu-Olivares (2007)
#'   is implemented.
#'   }
#'   
#'  \item{\code{se}} {
#'   When \code{se} is set to \code{"boot"} or \code{"bootstrap"} standard 
#'   errors are computed using the pairs bootstrap (Freedman, 1984).
#'   The standard errors are based on the standard deviation of successful 
#'   bootstrap replications. The \code{z-value} and \code{P(>|z|)} assume 
#'   the ratio of the coefficient estimate to the bootstrap standard 
#'   deviation approximates a normal  distribution.  Note, the Sargan 
#'   test statistic is calculated from the original sample and is not 
#'   a bootstrap-based estimate. When \code{se} is set to \code{standard}
#'   and \code{mising} is \code{listwise} standard errors for the 
#'   MIIV-2SLS coefficients are calculated using analytic expressions.
#'   For equations with categorical endogenous variables, standard
#'   errors of the MIIV-2SLS coefficients are calculated using the ADD MORE.........
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
#' Bollen, K. A. (1996).	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 
#' Bollen, K. A. (2001).	Two-stage	Least	Squares	and	Latent	Variable	
#' Models: Simultaneous	Estimation	and	Robustness	to	Misspecifications.
#' In	R.	Cudeck,	S.	Du	Toit,	and	D.	Sorbom	(Eds.),	Structural	
#' Equation	Modeling:	Present	and	Future,	A	Festschrift	in	Honor	of	Karl	
#' Joreskog	(pp. 119-138).	Lincoln,	IL: Scientific	Software.
#' 	
#' Bollen, K. A., & Maydeu-Olivares, A. (2007). A Polychoric Instrumental 
#' Variable (PIV) Estimator for Structural Equation Models with Categorical 
#' Variables. \emph{Psychometrika}, 72(3), 309.
#' 
#' Freedman, D. (1984). On Bootstrapping Two-Stage Least-Squares Estimates 
#' in Stationary Linear Models. \emph{The Annals of Statistics}, 
#' 12(3), 827–842. 
#' 
#' Greene, W. H. (2000). Econometric analysis. Upper Saddle River, N.J: 
#' Prentice Hall.
#' 
#' Hayashi, F. (2000). Econometrics. Princeton, NJ: Princeton University 
#' Press
#'
#' @example example/bollen1989-miive1.R
#' @example example/bollen1989-miive2.R
#' @example example/bollen1989-miive3.R
#' 
#' @seealso \link{MIIVsem}{miivs}
#' 
#' @keywords MIIV-2SLS MIIV PIV 2sls tsls instrument SEM two-stage least-squares
#'  
#' @export
miive <- function(model = model, 
                  data = NULL,  
                  instruments = NULL,
                  sample.cov = NULL, 
                  sample.mean = NULL, 
                  sample.nobs = NULL, 
                  sample.cov.rescale = TRUE, 
                  estimator = "2SLS", 
                  se = "standard", 
                  bootstrap = 1000L, 
                  missing = "listwise",
                  est.only = FALSE, 
                  var.cov = FALSE, 
                  var.cov.estimator = "ML",
                  miiv.check = TRUE, 
                  ordered = NULL, 
                  wald = FALSE){
  
  #-------------------------------------------------------# 
  # A few basic sanity checks for user-supplied covariance 
  # matrices and mean vectors. Also, rescale cov.matrix if
  # cov.rescale set to TRUE.
  #-------------------------------------------------------# 
  if(is.null(data)){
    
    if (!is.null(ordered)){
      stop(paste(
        "miive: raw data required when declaring factor variables.")
      )
    }
    
    if (!is.vector(sample.mean)){
      stop(paste(
        "miive: sample.mean must be a vector.")
      )
    }
    
    if (is.null(names(sample.mean))){
      stop(paste("miive: sample.mean vector must have names."))
    }
    
    if (!all.equal( names(sample.mean),
                    colnames(sample.cov), 
                    check.attributes = FALSE )){
      stop(paste(
        "miive: names of sample.mean vector",
        "and sample.cov matrix must match.")
      )
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
  # remove equations where there are not sufficient MIIVs
  #-------------------------------------------------------#
  underid   <- unlist(lapply(d, function(eq) {
    length(eq$MIIVs) < length(eq$IVobs)
  }))
  d.un <- d[underid]; d   <- d[!underid]
  
  #-------------------------------------------------------# 
  # Remove variables from data that are not in model 
  # syntax and preserve the original column ordering.
  #-------------------------------------------------------# 
  if(!is.null(data)){
    
    obs.vars <- unique(unlist(lapply(d,"[", c("DVobs", "IVobs", "MIIVs"))))
    
    if (any(!obs.vars %in% colnames(data))){
     stop(paste(
       "miive: model syntax contains variables not in data.")
     )
    }
    
    data <- data[,colnames(data) %in% obs.vars]
    data <- as.data.frame(data)
    
    # convert any variables listed in 
    # ordered to categorical
    data[,ordered]  <- lapply(
      data[,ordered, drop = FALSE], ordered
    )
  }
  
  #-------------------------------------------------------# 
  # If missing == "listwise"
  #-------------------------------------------------------# 
  if (!is.null(data) & missing == "listwise"){
    
    data <- data[stats::complete.cases(data),]
    
  } 

  #-------------------------------------------------------# 
  # Process data. See documentation of processRawData. 
  #-------------------------------------------------------# 
  g <- processData(data, 
                   sample.cov, 
                   sample.mean, 
                   sample.nobs, 
                   ordered, 
                   missing,
                   se,
                   pt )
  
  #-------------------------------------------------------# 
  # Add some fields to d and check for any problematic
  # cases prior to estimation.  
  #-------------------------------------------------------# 
  d <- lapply(d, function (eq) {
    
    eq$missing <- ifelse(
      any(g$var.nobs[c(eq$DVobs, eq$IVobs, eq$MIIVs)] < g$sample.nobs), 
      TRUE, 
      FALSE
    )
    
    eq$categorical <- ifelse(
      any(g$var.categorical[c(eq$DVobs, eq$IVobs, eq$MIIVs)]), 
      TRUE, 
      FALSE
    )
    
    eq$exogenous   <- ifelse(
      any(g$var.exogenous[c(eq$DVobs, eq$IVobs, eq$MIIVs)]), 
      TRUE, 
      FALSE
    )
    
    # For now throw an error if a categorical equation
    # contains an exogenous variable. 
    if (eq$categorical & eq$exogenous){
      
      stop(paste("miive: exogenous variables in",
                 "equations containing categorical",
                 "variables (including MIIVs)",
                 "are not currently supported."))
    }
    # if (eq$missing & eq$restricted){
    #   stop(paste("miive: Restrictions on coefficients in",
    #              "equations with missing values",
    #              "(including on the MIIVs)",
    #              "is not currently supported."))
    # }
    eq
  })
  
  #-------------------------------------------------------#  
  # Build Restriction Matrices.
  #-------------------------------------------------------#  
  r <- buildRestrictMat(d)
  
  d <- lapply(d, function (eq) {

    eq$restricted  <- ifelse(
      r$eq.restricted[eq$DVobs],
      TRUE,
      FALSE
    )
    # throw an error if a categorical equation 
    # includes restrictions
    if ((eq$categorical & eq$restricted) & se == "standard"){
      stop(paste("miive: When SE = 'standard', restrictions",
                 "on coefficients in equations containing",
                 "categorical variables (including MIIVs)",
                 "are not allowed. Remove restrictions or",
                 "use se = 'bootstrap'."))
    }
    eq
  })
  
  #-------------------------------------------------------#
  # Generate results object.
  #-------------------------------------------------------#
  results <- switch(
    estimator,
      "2SLS" = miive.2sls(d, g, r, est.only, se, missing, var.cov),
      # "GMM"  = miive.gmm(d, d.un, g, r, est.only, se), # Not implemented
      # In other cases, raise an error
      stop(
        paste(
          "Invalid estimator:", 
          estimator, 
          "Valid estimators are: 2SLS"
        )
      )
  )
  
  #-------------------------------------------------------#
  # Point estimates for variance and covariance params
  #-------------------------------------------------------#
  if (var.cov){
    
    if(length(d.un) > 0){
      stop(paste("MIIVsem: variance covariance esimtation not",
                 "allowed in the presence of underidentified",
                 "equations."))
    }
    
    v <- list()
    v$coefficients <- estVarCovarCoefs(data, g, results$eqn, pt, ordered)
    
  } else {
    
    v <- NULL
    
  }
  

  #-------------------------------------------------------#
  # Boostrap and substitute closed form SEs with boot SEs
  #-------------------------------------------------------#
  
  if(se == "boot" | se == "bootstrap"){
    
    # Do this outside of the bootstrap loop

    boot.results <- boot::boot(data, function(origData, indices){
      
      bsample <- origData[indices,]
      
      g.b <- processData(
        data = bsample, 
        sample.cov = NULL, 
        sample.mean = NULL, 
        sample.nobs = NULL, 
        ordered = ordered, 
        missing = missing,
        se = se,
        pt = pt
      )
      
      
      brep <- switch(
        estimator,
        "2SLS" = miive.2sls(d, g.b, r, est.only = TRUE, se, missing, var.cov),
        #"GMM"  = miive.gmm(d, g, r, est.only = TRUE), 
        # Not implemented
        # In other cases, raise an error
        stop(paste(
          "Invalid estimator:", 
          estimator, 
          "Valid estimators are: 2SLS")
        )
      )
      
      if (var.cov){
        
        v.b <- estVarCovarCoefs(bsample, g.b, brep$eqn, pt, ordered)
        
        c(brep$coefficients, v.b)
        
      } else {
        
        brep$coefficients
        
      }
      
    }, bootstrap)
    
    # Replace the estimated variances of estimates 
    # with the boostrap estimates
    boot.mat <- boot.results$t[
      stats::complete.cases(boot.results$t),
    ]
    
    results$bootstrap.true <- nrow(boot.mat)
    
    # Need to determine which columns are constants
    if (results$bootstrap.true > 0){
      
      nearZero <- apply(boot.mat , 2, stats::var, na.rm=TRUE) < 1e-16
      
    } else {
      
      nearZero <- FALSE
      
    }
    
    coefCov <- stats::cov(boot.mat)
    
    if(any(nearZero)) {
      coefCov[nearZero,nearZero] <- round(
        coefCov[nearZero,nearZero], 1e-10
      )
    }
    
    rownames(coefCov) <- colnames(coefCov) <- c(
      names(results$coefficients), 
      if (var.cov) names(v$coefficients) else NULL
    )
    
    results$eqn <- lapply(results$eqn, function(eq) {
      eq$coefCov <- coefCov[
          names(eq$coefficients),
          names(eq$coefficients),
          drop = FALSE
      ]
      eq
    })
    
    results$coefCov <- coefCov
    # Store the boot object as a part of the result object. 
    # This is useful for calculating CIs or
    # other bootstrap postestimation.
    results$boot <- boot.results

  } # end bootstrap
  

  if (missing == "twostage" & se == "standard" & var.cov == TRUE){
    
    coefCov <- estTwoStageML(g,v,results$eqn,pt)
    
    results$eqn <- lapply(results$eqn, function(eq) {
      eq$coefCov <- coefCov[
          names(eq$coefficients),
          names(eq$coefficients),
          drop = FALSE
      ]
      eq
    })
    
    results$coefCov <- coefCov
    v$coefcov <- coefCov[
         names(v$coefficients), 
         names(v$coefficients), 
         drop = FALSE
    ]
  }
  
  # assemble return object
  results$model          <- model
  results$estimator      <- estimator
  results$se             <- se
  results$missing        <- missing
  results$bootstrap      <- bootstrap
  results$call           <- match.call()
  results$ordered        <- ordered
  results$wald           <- wald
  results$eqn.unid       <- d.un
  results$r              <- r
  results$v              <- v

  
  class(results)  <- "miive"
  
  results
}

