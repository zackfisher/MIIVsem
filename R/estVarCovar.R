#' Estimate the variance and covariance parameters
#' 
#' @param data
#' @param vcov.model
#' @param ordered
#' @param var.cov.estimator
#'@keywords internal
estVarCovar <- function(data, vcov.model, ordered,var.cov.estimator){
  
  v <- list()
  
  v$coefficients  <- tryCatch(
    {
      pe <- lavaan::parameterEstimates(
        lavaan::sem(
          vcov.model,
          estimator = var.cov.estimator,
          se = "none",
          data,
          ordered = ordered
        )
      )
      pe <- pe[pe$op == "~~", , drop = FALSE]
      v.coefficients        <- pe[,"est"]
      names(v.coefficients) <- paste0(pe$lhs,"~~",pe$rhs)
      v.coefficients
    },
    error = function(cond) 
      { 
        zt <- lavaan::lavaanify(vcov.model)
        zt <- zt[zt$op == "~~", , drop = FALSE]
        v.coefficients        <- rep(NA, nrow(zt))
        names(v.coefficients) <- paste0(zt$lhs,"~~",zt$rhs)
        v.coefficients
      },
    warning = function(cond) 
      { 
        zt <- lavaan::lavaanify(vcov.model)
        zt <- zt[zt$op == "~~", , drop = FALSE]
        v.coefficients        <- rep(NA, nrow(zt))
        names(v.coefficients) <- paste0(zt$lhs,"~~",zt$rhs)
        v.coefficients
      },
    finally={}
  )    
  
  return(v)
}