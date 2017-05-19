#' return a dataframe of parameter estimates for a fitted model.
#' 
#' @param x An object of class miive
#' @param v A list containing variance snd covariance parameter information.
#' @param sarg Logical. Should Sargan test results be included.
#'  
#' @export
estimatesTable <- function(x, v = NULL, sarg = FALSE){
  
  # measurement equations
  meas.eqns <-  unlist(lapply(x$eqn, function(eq){
    ifelse(eq$EQmod == "measurement", TRUE, FALSE)
  }))
  
  # structural equations
  str.eqns <-  unlist(lapply(x$eqn, function(eq){
    ifelse(eq$EQmod == "regression", TRUE, FALSE)
  }))
  
  if (!is.null(x$boot)){
    
    zeroVariance <- function(x, tol = .Machine$double.eps ^ 0.5) {
      if (length(x) == 1) return(TRUE)
      x <- range(x) / mean(x)
      isTRUE(all.equal(x[1], x[2], tolerance = tol))
    }
    
    getAllCI <- function(z, type = type, w) {
      
      b.ci <- rbind(utils::tail(unlist(
        boot::boot.ci(z, index = w, type = type)[-(1:3)]
      ),2))
      
      rownames(b.ci) <- names(boot::boot.ci(z, index = w, type = type)$t0)
      colnames(b.ci) <- c("lwr","upr")
      
      b.ci
    }
    
    
    b.ci <- do.call("rbind",
      lapply(1:ncol(x$boot$t),function(col){
        if (zeroVariance(x$boot$t[,col])){
          tmp <- matrix(0,1,2)
          colnames(tmp) <- c("lwr","upr")
          rownames(tmp) <- names(x$boot$t0)[col]
          tmp
        } else {
          getAllCI(z = x$boot, type = x$boot.ci, w = col)
        }
      })
    )
    
  }
  
  
  
  
  # vector of named SEs
  se.diag <- sqrt(diag(x$coefCov))
  
  meas.coef.mat <- data.frame(
    "lhs" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      # if(eq$categorical){
      #   rep(eq$IVlat, length(eq$coefficients))
      # } else {
      #   c(eq$DVlat, rep(eq$IVlat, length(eq$coefficients)-1))
      # }
      if(eq$categorical){
        eq$IVlat
      } else {
        c(eq$DVlat, eq$IVlat)
      }
    })),
    "op"  = unlist(lapply(x$eqn[meas.eqns], function(eq){
      if(eq$categorical){
        rep("=~", length(eq$coefficients))
      } else {
        c("~1", rep("=~", length(eq$coefficients)-1))
      }
    })),
    "rhs" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      # if(eq$categorical){
      #   eq$DVlat
      # } else {
      #   c("", eq$DVlat)
      # }
      if(eq$categorical){
        rep(eq$DVlat, length(eq$IVlat))
      } else {
        c("", rep(eq$DVlat, length(eq$IVlat)))
      }
    })),
    "est" = unlist(lapply(x$eqn[meas.eqns], "[[", "coefficients")),
    "se" = unlist(lapply(x$eqn[meas.eqns], function(eq){
        se.diag[names(eq$coefficients)]
      })
    ),
    "z" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      if (is.null(x$boot)){
        eq$coefficients/se.diag[names(eq$coefficients)]
      } else {
         b.ci[names(eq$coefficients), 1]
      }
      })
    ),
    "pvalue" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      if (is.null(x$boot)){
        2*(stats::pnorm(abs(
          eq$coefficients/se.diag[names(eq$coefficients)]
        ), lower.tail=FALSE))
      } else {
        b.ci[names(eq$coefficients), 2]
      }
      })
    ),
    "sarg" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      rep(eq$sargan, length(eq$coefficients))})
    ),
    "sarg.df" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      rep(eq$sargan.df, length(eq$coefficients))})
    ),
    "sarg.p" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      rep(eq$sargan.p, length(eq$coefficients))})
    ), 
    # "eq" = unlist(lapply(x$eqn[meas.eqns], function(eq){
    #   rep(eq$EQnum, length(eq$coefficients))})
    # ), 
    stringsAsFactors = FALSE
  )
  
  if  (length(meas.eqns) > 0){
  # add a scaling indicator for each latent variable
    lv.si <- unique(do.call("rbind",
      lapply(x$eqn[meas.eqns], function(eq){
        cbind(eq$IVlat, eq$IVobs, eq$categorical)
      })
    ))
    


    
    # if there is a duplicated scaling indicator
    # in column 2 of lv.si, it should be replcad
    # by it's latent variable counterpart. 
    latVars <- unique(x$pt$lhs[x$pt$op=="=~"])
    higherOrderFactors <-  unique(x$pt$lhs[
      x$pt$op == "=~" &  x$pt$rhs %in% latVars
    ])
    
    if (length(higherOrderFactors) > 0){
      for ( i in 1:nrow(lv.si)){
        if (lv.si[i,1] %in% higherOrderFactors){
          lv.si[i,2]<- x$pt$rhs[x$pt$op == "=~" & x$pt$lhs == lv.si[i,1]][1]
        }
      }
    }

    
    if (!is.null(lv.si)){
      meas.coef.mat <- 
        rbind(do.call("rbind", apply(lv.si, 1, function(x){
          if (x[3]){
            data.frame(
              "lhs"     = x[1],
              "op"      = "=~",
              "rhs"     = x[2],
              "est"     = 1,
              "se"      = 0,
              "z"       = NA,
              "pvalue"  = NA,
              "sarg"    = NA,
              "sarg.df" = NA,
              "sarg.p"  = NA, 
              #"eq"      = NA, 
              stringsAsFactors = FALSE
            )
          } else {
            data.frame(
              "lhs"     = c(x[1],x[2]),
              "op"      = c("=~","~1"),
              "rhs"     = c(x[2],""),
              "est"     = c(1,0),
              "se"      = c(0,0),
              "z"       = NA,
              "pvalue"  = NA,
              "sarg"    = NA,
              "sarg.df" = NA,
              "sarg.p"  = NA, 
              #"eq"      = NA, 
              stringsAsFactors = FALSE
            )
          }
        })), meas.coef.mat)
    }
  }
  
  if(all(!is.null(x$sample.cov) & x$sample.mean == 0)){
    meas.coef.mat[meas.coef.mat$op == "~1","se"] <- 0
      meas.coef.mat[
        meas.coef.mat$op == "~1", 
        c("se","z","pvalue","sarg", "sarg.p")
      ] <- NA
  }
  
  str.coef.mat <- data.frame(
    "lhs" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$DVlat, length(eq$coefficients))})
    ),
    "op"  = unlist(lapply(x$eqn[str.eqns], function(eq){
      if(eq$categorical){
        rep("~", length(eq$coefficients))
      } else {
        c("~1", rep("~", length(eq$coefficients)-1))
      }
    })),
    "rhs" = unlist(lapply(x$eqn[str.eqns], function(eq){
      if(eq$categorical){
        eq$IVlat
      } else {
        c("", eq$IVlat)
      }
    })),
    "est" = unlist(
      lapply(x$eqn[str.eqns], "[[", "coefficients")
    ),
    "se" = unlist(lapply(x$eqn[str.eqns], function(eq){
        se.diag[names(eq$coefficients)]
      })
    ),
    "z" = unlist(lapply(x$eqn[str.eqns], function(eq){
      if(is.null(x$boot)){
        eq$coefficients/se.diag[names(eq$coefficients)]
      } else {
         b.ci[names(eq$coefficients), 1]
      }
      })
    ),
    "pvalue" = unlist(lapply(x$eqn[str.eqns], function(eq){
      if(is.null(x$boot)){
        2*(stats::pnorm(abs(
          eq$coefficients/se.diag[names(eq$coefficients)]
        ), lower.tail=FALSE))
      } else {
        b.ci[names(eq$coefficients), 2]
      }
      })
    ),
    "sarg" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$sargan, length(eq$coefficients))})
    ),
    "sarg.df" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$sargan.df, length(eq$coefficients))})
    ),
    "sarg.p" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$sargan.p, length(eq$coefficients))})
    ), 
    # "eq" = unlist(lapply(x$eqn[str.eqns], function(eq){
    #   rep(eq$EQnum, length(eq$coefficients))})
    # ), 
    stringsAsFactors = FALSE
  )
  
  if (!is.null(x$v)){
    v <- x$v
    vcov.names   <- names(v$coefficients)
    vcov.coefCov <- v$coefCov[
      vcov.names,
      vcov.names
    ]

    vcov.coef.mat <- data.frame(
      "lhs" = do.call(rbind, strsplit(vcov.names, "~~"))[,1],
      "op"  = "~~",
      "rhs" = do.call(rbind, strsplit(vcov.names, "~~"))[,2],
      "est" = v$coefficients,
      "se"  = if(is.null(v$coefCov)) NA else 
        sqrt(diag(vcov.coefCov)),
      "z" = if(is.null(v$coefCov)) {
        NA
       } else if (is.null(x$boot)){
         v$coefficients/sqrt(diag(vcov.coefCov))
       } else {
         b.ci[names(v$coefficients), 1]
       },
      "pvalue" = if(is.null(v$coefCov)) {
        NA
       } else if (is.null(x$boot)){
         2*(stats::pnorm(abs(
         v$coefficients/sqrt(diag(vcov.coefCov))), 
         lower.tail=FALSE))
       } else {
         b.ci[names(v$coefficients), 2]
       },
      "sarg" = NA,
      "sarg.df" = NA, 
      "sarg.p" = NA,
      #"eq" = NA,
      stringsAsFactors = FALSE
    )

  } else {
    
    vcov.coef.mat <- NULL
    
  }
  
  parTab <- rbind(
    str.coef.mat,
    meas.coef.mat,
    vcov.coef.mat
  )
  
  parTab[is.infinite(parTab[,"z"]),c("se","z","pvalue")] <- NA
  
  parTab <- parTab[order(parTab$op, parTab$lhs),] 
  
  if(!is.null(x$boot)){
    colnames(parTab)[colnames(parTab) ==      "z"] <- "lower"
    colnames(parTab)[colnames(parTab) == "pvalue"] <- "upper"
  }
  
  rownames(parTab) <- NULL
  
  if (!sarg){
    parTab$sarg <- parTab$sarg.df <- parTab$sarg.p <- NULL
  }

  
  return(parTab)
}