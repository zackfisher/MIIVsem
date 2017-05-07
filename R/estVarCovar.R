#' estimate the variance and covariance parameters
#'@keywords internal
estVarCovar <- function(data, 
                        g, 
                        results, 
                        pt, 
                        ordered, 
                        se, 
                        missing, 
                        var.cov.estimator){
  
  
  # if (!is.null(x$v)){
  #   v <- x$v
  #   vcov.names   <- names(v$coefficients)
  #   vcov.coefCov <- x$coefCov[
  #     colnames(x$coefCov) %in% vcov.names,
  #     colnames(x$coefCov) %in% vcov.names
  #     ]
  #   
  #   vcov.coef.mat <- data.frame(
  #     "lhs" = do.call(rbind, strsplit(vcov.names, "~~"))[,1],
  #     "op"  = "~~",
  #     "rhs" = do.call(rbind, strsplit(vcov.names, "~~"))[,2],
  #     "est" = v$coefficients,
  #     "se"  = if(dim(vcov.coefCov)[1] == 0) NA else 
  #       sqrt(diag(vcov.coefCov)),
  #     "z" = if(dim(vcov.coefCov)[1] == 0) NA else 
  #       v$coefficients/sqrt(diag(vcov.coefCov)),
  #     "pvalue" = if(dim(vcov.coefCov)[1] == 0) NA else  
  #       2*(stats::pnorm(abs(
  #         v$coefficients/sqrt(diag(vcov.coefCov))), 
  #         lower.tail=FALSE)),
  #     "sarg" = NA,
  #     "sarg.df" = NA, 
  #     "sarg.p" = NA,
  #     "eq" = NA,
  #     stringsAsFactors = FALSE
  #   )

  v <- list()
  
  # if there are categorical variables
  # and the estimator was left at the
  # default of "ML" set to  "DWLS"
  if(!is.null(ordered) & var.cov.estimator == "ML"){
    var.cov.estimator <- "DWLS"
  }
  
  
  
  

  v$coefficients  <- tryCatch(
    {
      if(!is.null(data)){
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
      } else {
        pe <- lavaan::parameterEstimates(
          lavaan::sem(
            vcov.model,
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
        zt <- lavaan::lavaanify(vcov.model)
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
  
  return(v)
}