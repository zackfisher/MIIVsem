#' return a dataframe of parameter estimates for a 
#' fitted model.
#' 
#' @param x an object of class miive
#' 
#' @export
estimatesTable <- function(x, v = NULL){
  
  #if(class(x) == "miive") { e <- x$eqn; v <- x$v } else { e <- x; }
  
  # measurement equations
  meas.eqns <-  unlist(lapply(x$eqn, function(eq){
    ifelse(eq$EQmod == "measurement", TRUE, FALSE)
  }))
  
  # structural equations
  str.eqns <-  unlist(lapply(x$eqn, function(eq){
    ifelse(eq$EQmod == "structural", TRUE, FALSE)
  }))
    
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
      sqrt(diag(eq$coefCov))})
    ),
    "z" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      eq$coefficients/sqrt(diag(eq$coefCov))})
    ),
    "pvalue" = unlist(lapply(x$eqn[meas.eqns], function(eq){
      2*(stats::pnorm(abs(eq$coefficients/sqrt(diag(eq$coefCov))), lower.tail=FALSE))})
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
  
  # add a scaling indicator for each latent variable
  
  if(any(meas.eqns)){
    
    tmp <- unique(
      do.call("rbind",
              lapply(x$eqn[meas.eqns], function(eq){
                if(!eq$categorical){
                     c(eq$IVlat, eq$IVobs)
                }
              })
      )
    )
    
    if (!is.null(tmp)){
      meas.coef.mat <- 
        rbind(do.call("rbind", apply(tmp, 1, function(x){
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
        })), meas.coef.mat)
    }
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
      sqrt(diag(eq$coefCov))})
    ),
    "z" = unlist(lapply(x$eqn[str.eqns], function(eq){
      eq$coefficients/sqrt(diag(eq$coefCov))})
    ),
    "pvalue" = unlist(lapply(x$eqn[str.eqns], function(eq){
      2*(stats::pnorm(abs(eq$coefficients/sqrt(diag(eq$coefCov))), lower.tail=FALSE))})
    ),
    "sarg" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$sargan, length(eq$coefficients))})
      #c(eq$sargan, rep(NA, length(eq$coefficients)-1))})
    ),
    "sarg.df" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$sargan.df, length(eq$coefficients))})
      #c(eq$sargan.df, rep(NA, length(eq$coefficients)-1))})
    ),
    "sarg.p" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$sargan.p, length(eq$coefficients))})
      #c(eq$sargan.p, rep(NA, length(eq$coefficients)-1))})
    ), 
    "eq" = unlist(lapply(x$eqn[str.eqns], function(eq){
      rep(eq$EQnum, length(eq$coefficients))})
    ), 
    stringsAsFactors = FALSE
  )
  
  if (!is.null(x$v)){
    v <- x$v
    vcov.names   <- names(v$coefficients)
    vcov.coefCov <- x$coefCov[
      colnames(x$coefCov) %in% vcov.names,
      colnames(x$coefCov) %in% vcov.names
    ]

    vcov.coef.mat <- data.frame(
      "lhs" = do.call(rbind, strsplit(vcov.names, "~~"))[,1],
      "op"  = "~~",
      "rhs" = do.call(rbind, strsplit(vcov.names, "~~"))[,2],
      "est" = v$coefficients,
      "se"  = if(dim(vcov.coefCov)[1] == 0) NA else 
        sqrt(diag(vcov.coefCov)),
      "z" = if(dim(vcov.coefCov)[1] == 0) NA else 
        v$coefficients/sqrt(diag(vcov.coefCov)),
      "pvalue" = if(dim(vcov.coefCov)[1] == 0) NA else  
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