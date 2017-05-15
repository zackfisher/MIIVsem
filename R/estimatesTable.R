#' return a dataframe of parameter estimates for a fitted model.
#' 
#' @param x an object of class miive
#' @param v a list containing variance snd covariance parameter information.
#'  
#' @export
estimatesTable <- function(x, v = NULL){
  
  # measurement equations
  meas.eqns <-  unlist(lapply(x$eqn, function(eq){
    ifelse(eq$EQmod == "measurement", TRUE, FALSE)
  }))
  
  # structural equations
  str.eqns <-  unlist(lapply(x$eqn, function(eq){
    ifelse(eq$EQmod == "regression", TRUE, FALSE)
  }))
  
  # vector of named SEs
  se.diag <- sqrt(diag(x$coefCov))
  
  meas.coef.mat <- data.frame(
    "lhs" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      if(eq$categorical){
        rep(eq$IVlat, length(eq$coefficients))
      } else {
        c(eq$DVlat, rep(eq$IVlat, length(eq$coefficients)-1))
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
      if(eq$categorical){
        eq$DVlat
      } else {
        c("", eq$DVlat)
      }
    })),
    "est" = unlist(lapply(x$eqn[meas.eqns], "[[", "coefficients")),
    "se" = unlist(lapply(x$eqn[meas.eqns], function(eq){
        se.diag[names(eq$coefficients)]
      })
    ),
    "z" = unlist(lapply(x$eqn[meas.eqns], function(eq){
        eq$coefficients/se.diag[names(eq$coefficients)]
      })
    ),
    "pvalue" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      2*(stats::pnorm(abs(
        eq$coefficients/se.diag[names(eq$coefficients)]
      ), lower.tail=FALSE))})
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
    "eq" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      rep(eq$EQnum, length(eq$coefficients))})
    ), 
    stringsAsFactors = FALSE
  )
  
  if  (length(meas.eqns) > 0){
  # add a scaling indicator for each latent variable
    lv.si <- unique(do.call("rbind",
      lapply(x$eqn[meas.eqns], function(eq){
        c(eq$IVlat, eq$IVobs, eq$categorical)
      })
    ))
    
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
              "eq"      = NA, 
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
              "eq"      = NA, 
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
      eq$coefficients/se.diag[names(eq$coefficients)]
      })
    ),
    "pvalue" = unlist(lapply(x$eqn[str.eqns], function(eq){
      2*(stats::pnorm(abs(
        eq$coefficients/se.diag[names(eq$coefficients)]
      ), lower.tail=FALSE))})
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
    "eq" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$EQnum, length(eq$coefficients))})
    ), 
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
      "z" = if(is.null(v$coefCov)) NA else 
        v$coefficients/sqrt(diag(vcov.coefCov)),
      "pvalue" = if(is.null(v$coefCov)) NA else  
        2*(stats::pnorm(abs(
          v$coefficients/sqrt(diag(vcov.coefCov))), 
          lower.tail=FALSE)),
      "sarg" = NA,
      "sarg.df" = NA, 
      "sarg.p" = NA,
      "eq" = NA,
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
  
  rownames(parTab) <- NULL
  
  return(parTab)
}