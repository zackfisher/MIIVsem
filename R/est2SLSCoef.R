#' estimate the 2SLS point estimates
#'@keywords internal
est2SLSCoef <- function(d = d, r = NULL, poly.mat = NULL, sscp.mat = NULL){
  
  # Construct the following matrices:
  # XY1: A vector of crossproducts of fitted values from the first stage
  #      regression and the dependent variables. The crossproducts
  #      for all regressions are stacked into one vector
  # XX1: A block diagnonal matrix of of crossproducts of fitted values
  #      from the first stage regressions for all equations.
  # ZV:  A block diagonal matrix of crossproducts between the instuments
  #      and independent variables for all equations
  
   sVV <- lapply(d, function(eq){
    if(eq$categorical){
      matInv(
        poly.mat[
          eq[["MIIVs"]], eq[["MIIVs"]], drop = FALSE
        ]
      )
    } else {
      matInv(
        sscp.mat[
          c("1",eq[["MIIVs"]]), c("1",eq[["MIIVs"]]), drop = FALSE
        ]
      )
    }
  })
  
  ZV <- lapply(d, function(eq){
    if(eq$categorical){
      poly.mat[
        eq[["IVobs"]], eq[["MIIVs"]], drop = FALSE
      ]
    } else {
      sscp.mat[
        c("1",eq[["IVobs"]]), c("1",eq[["MIIVs"]]), drop = FALSE
      ]
    }
  })

  VY <- lapply(d, function(eq){
    if(eq$categorical){
      poly.mat[
        eq[["MIIVs"]], eq[["DVobs"]], drop = FALSE
      ]
    } else {
      sscp.mat[
        c("1",eq[["MIIVs"]]), eq[["DVobs"]], drop = FALSE
      ]
    }
  })
  
  XX1 <- mapply(function(ZV,sVV){
    
    ZV %*% sVV %*% t(ZV)
    
  },ZV,sVV, SIMPLIFY = FALSE)
  
  XY1 <- mapply(function(ZV,sVV,VY){
    
    ZV %*% sVV %*% VY
    
  },ZV,sVV,VY, SIMPLIFY = FALSE)
  
  
  ## coefficient names
  coef.names <- lapply(d, function(eq) { 
    paste0(eq$DVlat, "~",if(eq$categorical) eq$IVlat else c("1",eq$IVlat))
  })
  
  
  #------------------------------------------------------------------------#
  # MIIV-2SLS point estimates: unrestricted
  #------------------------------------------------------------------------#
  if (is.null(r$R)){ 
    
    
    d <- mapply(function (d, XX1, XY1, coef.names){
      
           d$coefficients <- as.numeric(solve(XX1, XY1))
           d$XX1          <- XX1
           names(d$coefficients) <-  coef.names; d
           
    }, d, XX1, XY1, coef.names, SIMPLIFY = FALSE)
    
    #coefficients <- unlist(lapply(d, "[[","coefficients"))
  
  #------------------------------------------------------------------------#
  # MIIV-2SLS point estimates: restricted
  #------------------------------------------------------------------------#
  } else { 
    
    R0   <- matrix(0, nrow(r$R), nrow(r$R))
    
    XX1F <- lavaan::lav_matrix_bdiag(XX1)
    
    M1   <- rbind(cbind(XX1F, t(r$R) ), 
                  cbind(r$R,     R0  ))
    
    M2   <- rbind(cbind(unlist(XY1)), 
                  r$q              )
    
    coefficients <- as.numeric(solve(M1,M2)[1:ncol(r$R),])
    
    # Add coefficients to equations list.
    coefList <- split(
      coefficients, 
      unlist(lapply(seq_along(d), function(x){
        if (d[[x]]$categorical){
          numCoefs <- length(d[[x]]$IVobs)
        } else {
          numCoefs <- length(d[[x]]$IVobs) + 1
        }
       rep(x,numCoefs)
      }))
    )
    
    names(coefList) <- rep("coefficients", length(coefList))
    names(XX1) <- rep("XX1", length(XX1))
    
    d  <- lapply(seq_along(d), function(x){
      u <- append(d[[x]],coefList[x])
      u <- append(u, XX1[x])
      names(u$coefficients) <- unlist(coef.names[x])
      u
    }) 
    
  }
  
  return(d)
}