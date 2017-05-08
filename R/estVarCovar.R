#' estimate the variance and covariance parameters
#'@keywords internal
estVarCovar <- function(data,
                        g,
                        eqns, 
                        pt, 
                        ordered, 
                        se, 
                        missing){
  
  var.cov.estimator <- "ML"
  var.cov.missing <- "FIML"
  
  v <- list()
  
  # if there are categorical variables
  # and the estimator was left at the
  # default of "ML" set to  "DWLS"
  if(!is.null(ordered) & var.cov.estimator == "ML"){
    var.cov.estimator <- "DWLS"
    var.cov.missing <- "listwise"
  }
  
  # fill the parTable 
  pt <- fillParTable(eqns, pt)
  
  # generate the model for variance covariance point 
  # estimates
  var.cov.model <- buildVarCovSyntax(pt)
  
  # estimate the variance covariance point estimates
  v$coefficients  <- tryCatch(
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

  #-------------------------------------------------------#
  # missing 
  #-------------------------------------------------------#
  if (se != "boot" & se != "bootstrap"){
    
    if (missing == "savalei"){
      
      pt <- fillParTable(eqns, pt, v)
      
      var.cov.se.model <- buildMissingSyntax(pt)
      
      # FIXME: remove suppress warnings.
      # temporarily added suppress warnings here as 
      # lavaan throws warnings re convergence 
      # when iter.max is set to "none"
      suppressWarnings( 
        var.cov.se.fit <- lavaan::lavaan(
          model = var.cov.se.model,
          sample.cov = g$sample.cov, 
          sample.mean = g$sample.mean,
          sample.nobs = g$sample.nobs,
          control=list(iter.max="none")
        )
      )
      
      DA <- unclass(lavaan::lavInspect(var.cov.se.fit, "delta"))
      DA <- DA[c(grep("~~", rownames(DA)),grep("~1", rownames(DA))), ]
      
      mi.cov   <- inspect(var.cov.se.fit, "cov.ov")
      mi.cov   <- mi.cov[rownames(g$sample.cov), colnames(g$sample.cov)]
      mi.mean  <- inspect(var.cov.se.fit, "mean.ov")
      mi.mean  <- mi.mean[names(g$sample.mean)]
      
      D        <- buildDuplication(length(g$sample.mean))
      sigInv   <- solve(mi.cov)
      meanDiff <- g$sample.mean - mi.mean
      
      UL <- t(D) %*% 
        kronecker(
          sigInv, 
          sigInv %*% (g$sample.cov + tcrossprod(meanDiff)) %*% sigInv - .5*sigInv
        ) %*% D
      
      LL <- kronecker(
        sigInv, 
        t(meanDiff) %*% sigInv 
      ) %*% D
      
      JB <- cbind(rbind(UL, LL), rbind(t(LL), sigInv))
      
      ## Savalei (2010, Equation 7)
      coefCov <- solve(t(DA) %*% JB %*% DA) %*% 
        t(DA) %*% JB %*% g$asymptotic.cov.sat %*% JB %*% 
        DA %*% solve(t(DA) %*% JB %*% DA)
      
      colnames(coefCov) <- rownames(coefCov) <- colnames(DA)
      
      v$coefCov <- coefCov
      
    } else {
      
      v$coefCov <- NULL
      
    }
    
  }

  
  return(v)
}