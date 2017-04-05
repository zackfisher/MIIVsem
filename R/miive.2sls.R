#' Two-stage least square estimator for a system of equations
#' 
#' @param d a list containing equation information
#' @param g a list containing data information
#' @param r a list containing coefficient restrictions 
#' @param est.only should we only calculate coefficient estimates
#' 
#' @details 
#' MIIV estimation using estimation functions. An 
#' estimation function retuns a list containing the
#' following elements:
#'
#' coefficients - a vector of estimated coefficients
#' coefCov      - variance covariance matrix of the estimates
#'                (optional, depending on the est.only argument)
#' residCov     - variance covariance matrix of the equation
#'                disturbances (optional, depending on the 
#'                est.only argument).
#' eqn -          a list of miiv equations and equation level
#'                estimation results. fields added include:
#'                  coefficients: coefficients
#'                  sigma       : equation sigma^2
#'                  sargan      : Sargan's test statistic
#'                  sargan.df   : df freedom for Sargan's
#'
#'@keywords internal
miive.2sls <- function(d, g, r, est.only){
    
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
      chol2inv(chol(g$sample.polychoric[eq[["MIIVs"]], eq[["MIIVs"]], drop = FALSE]))
    } else {
      chol2inv(chol(g$sample.sscp[c("1",eq[["MIIVs"]]), c("1",eq[["MIIVs"]]), drop = FALSE]))
    }
  })
  
  ZV <- lapply(d, function(eq){
    if(eq$categorical){
      g$sample.polychoric[eq[["IVobs"]], eq[["MIIVs"]], drop = FALSE]
    } else {
      g$sample.sscp[c("1",eq[["IVobs"]]), c("1",eq[["MIIVs"]]), drop = FALSE]
    }
  })

  VY <- lapply(d, function(eq){
    if(eq$categorical){
      g$sample.polychoric[eq[["MIIVs"]], eq[["DVobs"]], drop = FALSE]
    } else {
      g$sample.sscp[c("1",eq[["MIIVs"]]), eq[["DVobs"]], drop = FALSE]
    }
  })
  
  XX1 <- mapply(function(ZV,sVV){ZV %*% sVV %*% t(ZV)},ZV,sVV)
  XY1 <- mapply(function(ZV,sVV,VY){ZV %*% sVV %*% VY},ZV,sVV,VY)
  
  if (is.null(r$R)){ 
    
    d <- mapply(function (d, XX1, XY1){
      d$coefficients <- as.numeric(solve(XX1, XY1))
      names(d$coefficients) <-  paste0(
        d$DVlat,"~", if(d$categorical) d$IVlat else c("1",d$IVlat)
      )
      d
    }, d, XX1, XY1, SIMPLIFY = FALSE)
    
    coefficients <- unlist(lapply(d, "[[","coefficients"))
  
    
  } else { 
    
    R0 <- matrix(0, nrow(r$R), nrow(r$R))
    
    coefficients <- as.numeric(
      (solve(rbind(cbind(lavaan::lav_matrix_bdiag(XX1), t(r$R)), 
       cbind(r$R, R0))) %*% rbind(cbind(unlist(XY1)), r$q))[1:ncol(r$R),]
    )
    
    names(coefficients) <- unlist(lapply(d, function(eq) {
      paste0(eq$DVlat,"~", if(eq$categorical) eq$IVlat else c("1",eq$IVlat))
    }))
    
    # Add coefficients to equations list.
    coefList <- split(coefficients, unlist(lapply(seq_along(d), function(x){
      rep(x,(length(d[[x]]$IVobs) + ifelse(d[[x]]$categorical, 0, 1)))
    })))
    
    names(coefList) <- rep("coefficients", length(coefList))
    d  <- lapply(seq_along(d), function(x){append(d[[x]],coefList[x])}) 
    
  } 
  
  residCov <- NULL
  coefCov  <- NULL
  
  if(!est.only){

    if (!is.null(g$sample.cov)){

      d <- lapply(d, function(eq) {
        if (eq$categorical){
          eq$sigma <- NA
        } else {
          eq$sigma <- (g$sample.cov[eq$DVobs, eq$DVobs] + (t(eq$coefficients[-1]) %*%
                       g$sample.cov[c(eq$IVobs), c(eq$IVobs)] %*% eq$coefficients[-1]) -
                      (2 * g$sample.cov[eq$DVobs, c(eq$IVobs)] %*% eq$coefficients[-1]))
        }
        eq
      })

      if (is.null(r$R)){

        d <- mapply(function (d, XX1){
          if(is.na(d$sigma)){
            d$coefCov <- NA
          } else {
            d$coefCov <-solve(XX1 %*% t(solve(diag(rep(d$sigma, length(d$coefficients))))))
            rownames(d$coefCov) <- colnames(d$coefCov) <- names(d$coefficients)
          }
        d
        }, d, XX1, SIMPLIFY = FALSE)
        
        coefCov <- lavaan::lav_matrix_bdiag(lapply(d,"[[","coefCov"))
        rownames(coefCov) <- colnames(coefCov) <- names(coefficients)
        
      } else {
        
        sig  <- diag(unlist(lapply(d, function(eq)rep(eq$sigma, length(eq$coefficients)))))
        XX1F <- lavaan::lav_matrix_bdiag(XX1)
        coefCov <- solve(rbind(cbind((XX1F %*% t(solve(sig))), t(r$R)),
                               cbind(r$R, R0)))[1:nrow(XX1F), 1:nrow(XX1F)]
        rownames(coefCov) <- colnames(coefCov) <- names(coefficients)
      }
    } 
    
    if (!is.null(g$sample.polychoric)){
      
      d <- lapply(d, function(eq){
        if(eq$categorical) { # eq <- d[[1]]; mat <- g$sample.polychoric
          K <- buildCategoricalK(eq, g$sample.polychoric)
          eq$coefCov <-  K %*% g$asymptotic.cov %*% t(K)
          rownames(eq$coefCov) <- colnames(eq$coefCov) <- names(eq$coefficients)
        }
        eq
      })
      
    }
    

    # Calculate the residual covariance matrix for the full system of
    # equations based on covariance matrix input.

    # dvs <- unlist(lapply(d, "[[", "DVobs"))
    # B   <- diag(length(dvs))
    # colnames(B) <- rownames(B) <- dvs
    # idx <- do.call("rbind",lapply(d, function(eq){
    #   cbind(eq$DVobs, eq$IVobs, if (any(eq$categorical)) eq$coefficients else eq$coefficients[-1] )
    # }))
    # idx <- idx[idx[,2] %in% dvs, ,drop = FALSE]
    # B[idx[,2:1, drop = FALSE]] <- -1*as.numeric(idx[,3])
    # 
    # if (!is.null(g$sample.cov)){
    #   residCov <- t(B) %*% g$sample.cov[dvs,dvs] %*% B
    # }
    # 
    # rownames(residCov) <- colnames(residCov) <- unlist(lapply(d,"[","DVlat"))

    # Sargan's test (Hayashi, p. 228)
    # what about restrictions?
    d <- lapply(d, function(eq) {
      if (eq$categorical){
        eq$sargan    <- NA
        eq$sargan.df <- NA
        eq$sargan.p  <- NA
      } else {
        eq$sargan <-
          (t(g$sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] -
               g$sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*%
               eq$coefficients[-1]) %*%
             solve(g$sample.cov[eq$MIIVs,eq$MIIVs]) %*%
             (g$sample.cov[eq$MIIVs,eq$DVobs, drop = FALSE] -
                g$sample.cov[eq$MIIVs,eq$IVobs, drop = FALSE] %*%
                eq$coefficients[-1]) /  eq$sigma)*g$sample.nobs
        eq$sargan.df <- length(eq$MIIVs) - length(eq$IVobs)
        eq$sargan.p <- pchisq(eq$sargan, eq$sargan.df, lower.tail = FALSE)
        
      }
      eq
     })
    
  }
  g$coefficients <- coefficients
  g$coefCov      <- coefCov
  g$eqn          <- d
  return(g)
}

