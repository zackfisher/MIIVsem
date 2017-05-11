#' two-stage least square estimator for a system of equations
#' 
#' @param d a list containing equation information
#' @param g a list containing data matrices and characteristics
#' @param r a list containing coefficient restrictions 
#' @param est.only should we only calculate coefficient estimates
#' @param se se estimation
#' 
#'@keywords internal
miive.2sls <- function(d, g, r, est.only, se, missing, var.cov){
    
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
      chol2inv(chol(
        g$sample.polychoric[
          eq[["MIIVs"]], eq[["MIIVs"]], drop = FALSE
        ]
      ))
    } else {
      chol2inv(chol(
        g$sample.sscp[
          c("1",eq[["MIIVs"]]), c("1",eq[["MIIVs"]]), drop = FALSE
        ]
      ))
    }
  })
  
  ZV <- lapply(d, function(eq){
    if(eq$categorical){
      g$sample.polychoric[
        eq[["IVobs"]], eq[["MIIVs"]], drop = FALSE
      ]
    } else {
      g$sample.sscp[
        c("1",eq[["IVobs"]]), c("1",eq[["MIIVs"]]), drop = FALSE
      ]
    }
  })

  VY <- lapply(d, function(eq){
    if(eq$categorical){
      g$sample.polychoric[
        eq[["MIIVs"]], eq[["DVobs"]], drop = FALSE
      ]
    } else {
      g$sample.sscp[
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
           
           names(d$coefficients) <-  coef.names; d
           
    }, d, XX1, XY1, coef.names, SIMPLIFY = FALSE)
    
    coefficients <- unlist(lapply(d, "[[","coefficients"))
  
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
    names(coefficients) <- unlist(coef.names)
    
    # Add coefficients to equations list.
    coefList <- split(
      coefficients, 
      unlist(lapply(seq_along(d), function(x){
        rep(x,length(coef.names[[x]]))
      }))
    )
    
    names(coefList) <- rep("coefficients", length(coefList))
    
    d  <- lapply(seq_along(d), function(x){
      append(d[[x]],coefList[x])
    }) 
    
  } 
  
  residCov <- NULL
  coefCov  <- NULL
  
  #------------------------------------------------------------------------#
  # MIIV-2SLS se and overidentification tests
  # Not run during bootstrap resampling
  #------------------------------------------------------------------------#
  if(!est.only){ 

    #----------------------------------------------------------------------#
    # If there are any continous DV equations
    #----------------------------------------------------------------------#
    if (!is.null(g$sample.cov)){
      
        not.cat <- unlist(lapply(d, function(eq){
            rep(!eq$categorical, length(eq$coefficients))
        }))

      d <- lapply(d, function(eq) {
        if (eq$categorical){
          eq$sigma <- NA
        } else {
          eq$sigma <- 
            (g$sample.cov[eq$DVobs, eq$DVobs] + 
               (t(eq$coefficients[-1]) %*%
             g$sample.cov[c(eq$IVobs), c(eq$IVobs)] %*% 
               eq$coefficients[-1]) -
            (2 * g$sample.cov[eq$DVobs, c(eq$IVobs)] %*% 
               eq$coefficients[-1]))
        }
        eq
      })
    
      #--------------------------------------------------------------------#
      # missing == "listwise"
      #--------------------------------------------------------------------#
      if (missing == "listwise"){
        #------------------------------------------------------------------#
        # MIIV-2SLS SEs: unrestricted
        #------------------------------------------------------------------#
        if (is.null(r$R)){
          d <- mapply(function (d, XX1,coef.names){
            if(is.na(d$sigma)){
              d$coefCov <- NA
            } else {
              d$coefCov <-solve(
                XX1 %*% t(
                  solve(diag(rep(d$sigma, length(d$coefficients))))
                )
              )
              rownames(d$coefCov) <- colnames(d$coefCov) <- coef.names
            }
          d
          }, d, XX1,coef.names, SIMPLIFY = FALSE)
          
          coefCov <- lavaan::lav_matrix_bdiag(lapply(d,"[[","coefCov"))
          rownames(coefCov) <- colnames(coefCov) <- unlist(coef.names)
        
        #------------------------------------------------------------------#
        # MIIV-2SLS SEs: restricted
        #------------------------------------------------------------------#
        } else {
        
          sig  <- diag(unlist(lapply(d, function(eq){
              if(!eq$categorical){
                rep(eq$sigma, length(eq$coefficients))
              }
          })))
          
          
          XX1F <- lavaan::lav_matrix_bdiag(
            mapply(function (d, XX1){
              if(!d$categorical){
                XX1
              } else {
                matrix(0,0,0)
              }
            }, d, XX1, SIMPLIFY = FALSE)
          )
          
          R1 <- r$R[,not.cat,drop=FALSE]
          
          coefCov <- solve(
            rbind( cbind(XX1F %*% t(solve(sig)), t(R1)),
                   cbind(R1,                     R0   )
            )
          )[1:nrow(XX1F), 1:nrow(XX1F)]
          
          
          rownames(coefCov) <- colnames(coefCov) <- unlist(coef.names)[not.cat]
          
          d <- lapply(d, function(eq) {
            if (eq$categorical) {
              eq$coefCov <- NA
            } else {
              eq$coefCov <- coefCov[
                names(eq$coefficients),
                names(eq$coefficients),
                drop = FALSE
              ]
            }
            eq
          })
        }
        
      #--------------------------------------------------------------------#
      # missing == "twostage" and var.cov == FALSE
      #--------------------------------------------------------------------#        
      } else if (missing == "twostage" & var.cov == FALSE) {
        
        acov.names <- colnames(g$asymptotic.cov.sat)
        K <- buildDeltaK(d, g$sample.cov, g$sample.mean, acov.names)
        coefCov <- K %*% g$asymptotic.cov.sat %*% t(K)

        
        if (!is.null(r$R)){
          
          R1 <- r$R[,not.cat,drop=FALSE]
          
          coefCov <- rbind( 
            cbind(coefCov,                t(R1)),
            cbind(R1,                     R0   )
          )[1:nrow(coefCov), 1:nrow(coefCov)]
          
          
        }
          
        rownames(coefCov) <- colnames(coefCov) <- unlist(coef.names)[not.cat]
          
        d <- lapply(d, function(eq) {
          if (eq$categorical) {
            eq$coefCov <- NA
          } else {
            eq$coefCov <- coefCov[
              names(eq$coefficients),
              names(eq$coefficients),
              drop = FALSE
            ]
          }
          eq
        })
        
      }
        
    }
    
    if (!is.null(g$sample.polychoric)){
      
      if (se == "standard"){
        
        d <- lapply(d, function(eq){
          
          if(eq$categorical) { 
            
            K <- buildCategoricalK(eq, g$sample.polychoric)
            
            eq$coefCov <-  K %*% g$asymptotic.cov %*% t(K)
            
            rownames(eq$coefCov) <- colnames(eq$coefCov) <- 
            names(eq$coefficients)
            
          }
          eq
        })
        
        coefCov <- lavaan::lav_matrix_bdiag(lapply(d,"[[","coefCov"))
        rownames(coefCov) <- colnames(coefCov) <- unlist(coef.names)
      }
      
    }
    
    # Sargan's test (Hayashi, p. 228)
    d <- lapply(d, function(eq) {
      
      if (eq$categorical | (length(eq$MIIVs) - length(eq$IVobs)) < 1 ){
        
        eq$sargan <- NA; eq$sargan.df <- NA; eq$sargan.p  <- NA
        
      } else {
        
        eq$sargan.df <- length(eq$MIIVs) - length(eq$IVobs)
        
        eq$sargan <-
          ( 
            t( g$sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] -
               g$sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*%
               eq$coefficients[-1]) %*%
               solve(g$sample.cov[eq$MIIVs,eq$MIIVs] ) 
            %*%
             ( g$sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] -
               g$sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*%
               eq$coefficients[-1] ) 
            /  eq$sigma 
          ) * g$sample.nobs
          
          eq$sargan.p <- stats::pchisq(
            
            eq$sargan, eq$sargan.df, lower.tail = FALSE
            
          )
      }
      eq
     })
    
  }
  
  g$coefficients <- coefficients
  g$coefCov      <- coefCov
  g$eqn          <- d
  
  return(g)
  
}

