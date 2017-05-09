#' estimate the variance and covariance parameters
#'@keywords internal
estVarCovarCoefs <- function(data, g, eqns, pt, ordered){
  
  var.cov.estimator <- "ML"
  var.cov.missing <- "FIML"
  
  # if there are categorical variables
  # and the estimator was left at the
  # default of "ML" set to  "DWLS"
  if(!is.null(ordered) & var.cov.estimator == "ML"){
    var.cov.estimator <- "DWLS"
    var.cov.missing   <- "listwise"
  }
  
  # fill the parTable 
  pt <- fillParTable(eqns, pt)
  
  # generate the model for variance covariance point 
  # estimates
  var.cov.model <- buildVarCovSyntax(pt)
  
  # estimate the variance covariance point estimates
  coefficients  <- tryCatch(
    {
      if(!is.null(data)){
        pe <- lavaan::parameterEstimates(
          lavaan::lavaan(
            var.cov.model,
            estimator = var.cov.estimator,
            missing = var.cov.missing,
            se = "none",
            data,
            ordered = ordered
          )
        )
        pe <- pe[pe$op == "~~", , drop = FALSE]
        v.coefficients        <- pe[,"est"]
        names(v.coefficients) <- paste0(pe$lhs,"~~",pe$rhs)
        v.coefficients
      } else {
        pe <- lavaan::parameterEstimates(
          lavaan::lavaan(
            var.cov.model,
            estimator = var.cov.estimator,
            se = "none",
            sample.cov  = g$sample.cov,
            sample.mean = g$sample.mean,
            sample.nobs = g$sample.nobs,
            ordered = ordered
          )
        )
        pe <- pe[pe$op == "~~", , drop = FALSE]
        v.coefficients        <- pe[,"est"]
        names(v.coefficients) <- paste0(pe$lhs,"~~",pe$rhs)
        v.coefficients
      }
    },
    error = function(cond) 
    { 
      zt <- lavaan::lavaanify(var.cov.model)
      zt <- zt[zt$op == "~~", , drop = FALSE]
      v.coefficients        <- rep(NA, nrow(zt))
      names(v.coefficients) <- paste0(zt$lhs,"~~",zt$rhs)
      v.coefficients
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

  return(coefficients)
}