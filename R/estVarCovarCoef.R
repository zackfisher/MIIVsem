#' estimate the variance and covariance parameters
#'@keywords internal
estVarCovarCoef <- function(data, g, eqns, pt, ordered,
                            vc.coef.estimator = "ML",
                            vc.coef.missing = "FIML"){
  
  # if there are categorical variables
  # and the estimator was left at the
  # default of "ML" set to  "ULS"
  if(!is.null(ordered) & vc.coef.estimator == "ML"){
    vc.coef.estimator <- "ULS"
    vc.coef.missing   <- "listwise"
  }
  
  # fill the parTable 
  pt <- fillParTable(eqns, pt)
  
  # generate the model for variance covariance point 
  # estimates
  var.cov.model <- buildVarCovSyntax(pt)
  
  # estimate the variance covariance point estimates
  v  <- tryCatch(
    {
      if(!is.null(data)){
  
        fit <- lavaan::lavaan(
                 var.cov.model,
                 data = data,
                 estimator = vc.coef.estimator,
                 missing = vc.coef.missing,
                 se = "none",
                 ordered = ordered)
        
        rsq <- lavaan::inspect(fit, "rsquare")
        pe  <- lavaan::parameterEstimates(fit)
        pe <- pe[pe$op == "~~", , drop = FALSE]
        coefficients        <- pe[,"est"]
        names(coefficients) <- paste0(pe$lhs,"~~",pe$rhs)
        list(coefficients = coefficients, 
             coefCov = NULL,
             rsquare = rsq,
             model = var.cov.model,
             pt = pt)
      } else {
        fit <- lavaan::lavaan(
                var.cov.model,
                estimator = vc.coef.estimator,
                missing = vc.coef.missing,
                se = "none",
                sample.cov  = g$sample.cov,
                sample.mean = g$sample.mean,
                sample.nobs = g$sample.nobs,
                sample.cov.rescale = FALSE,
                ordered = ordered)
        rsq <- lavaan::inspect(fit, "rsquare")
        pe <- lavaan::parameterEstimates(fit)
        pe <- pe[pe$op == "~~", , drop = FALSE]
        coefficients        <- pe[,"est"]
        names(coefficients) <- paste0(pe$lhs,"~~",pe$rhs)
        list(coefficients = coefficients,
             coefCov = NULL,
             rsquare = rsq,
             model = var.cov.model,
             pt = pt)
      }
    },
    error = function(cond) 
    { 
      zt <- lavaan::lavaanify(var.cov.model)
      zt <- zt[zt$op == "~~", , drop = FALSE]
      rsq <- rep(NA, length(eqns))
      coefficients        <- rep(NA, nrow(zt))
      names(coefficients) <- paste0(zt$lhs,"~~",zt$rhs)
      list(coefficients = coefficients, 
           coefCov = NULL,
           rsquare = rsq,
           model = var.cov.model,
           pt = pt)
    },
    # warning = function(cond) 
    #   { 
    #     zt <- lavaan::lavaanify(vcov.model)
    #     zt <- zt[zt$op == "~~", , drop = FALSE]
    #     v.coefficients        <- rep(NA, nrow(zt))
    #     names(v.coefficients) <- paste0(zt$lhs,"~~",zt$rhs)
    #     v.coefficients
    #   },
    finally={}
  )    

  return(v)
}