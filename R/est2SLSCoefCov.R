#' estimate the 2SLS coefficient covariance matrix
#'@keywords internal
est2SLSCoefCov <- function(d, 
                           poly.mat = NULL,
                           cov.mat  = NULL,
                           mean.vec = NULL, 
                           acov     = NULL, 
                           acov.sat = NULL,
                           r        = NULL){
  
  coef.names <- lapply(d, function(eq) { 
    paste0(eq$DVlat, "~",if(eq$categorical) eq$IVlat else c("1",eq$IVlat))
  })
  
  not.cat <- unlist(lapply(d, function(eq){
    rep(!eq$categorical, length(eq$coefficients))
  }))

  coefCov <- matrix(0, length(unlist(coef.names)),  
                       length(unlist(coef.names)) )
  
  rownames(coefCov) <- colnames(coefCov) <- unlist(coef.names)
  
  if (!is.null(acov.sat)) {
    
    # delta method standard errors

    acov.names <- colnames(acov.sat)

    K <- buildDeltaK(d, cov.mat, mean.vec, acov.names)

    coefCov <- K %*% acov.sat %*% t(K)

    rownames(coefCov) <- colnames(coefCov) <- unlist(coef.names)[not.cat]
    

  } else {
    
    
    d <- lapply(d, function(eq) {
    
      if (eq$categorical){
      
        eq$sigma   <- NA
        K <- buildCategoricalK(eq, poly.mat)
        eq$coefCov <-  K %*% acov %*% t(K)
        rownames(eq$coefCov) <- colnames(eq$coefCov) <- names(eq$coefficients)
      
      } else {
      
        eq$sigma <- 
          (cov.mat[eq$DVobs, eq$DVobs] + 
              (t(eq$coefficients[-1]) %*%
          cov.mat[c(eq$IVobs), c(eq$IVobs)] %*% 
            eq$coefficients[-1]) -
            (2 * cov.mat[eq$DVobs, c(eq$IVobs)] %*% 
            eq$coefficients[-1]))
      
        if (!is.null(r$R)){
          
          eq$coefCov <- NA
          
        } else {
          
          eq$coefCov <- solve(eq$XX1 %*% t(
            solve(diag(rep(eq$sigma, length(eq$coefficients))))))
            
          rownames(eq$coefCov) <- colnames(eq$coefCov) <- 
            paste0(eq$DVlat, "~", c("1",eq$IVlat))

        }
        
      }
      eq
    })
    
    
    
    if (is.null(r$R)){
      
      coefCov <- lavaan::lav_matrix_bdiag(lapply(d,"[[","coefCov"))
      rownames(coefCov) <- colnames(coefCov) <- unlist(coef.names)
      
    } else {
    
      SIG <- diag(unlist(lapply(d,function(eq){
        if (!eq$categorical) rep(eq$sigma,length(eq$coefficients))}))
      )
    
      XX1 <- lavaan::lav_matrix_bdiag(lapply(d, function (eq){
        if (!eq$categorical) eq$XX1 else matrix(0,0,0)
      }))
    
      R0 <- matrix(0, nrow(r$R), nrow(r$R))
      R1 <- r$R[,not.cat,drop=FALSE]
    
      coefCovR <- solve(
        rbind(cbind(XX1 %*% t(solve(SIG)), t(R1)),cbind(R1,R0))
      )[1:nrow(XX1), 1:nrow(XX1)]
    
      colnames(coefCovR) <- rownames(coefCovR) <- 
        unlist(coef.names)[not.cat]
      
      coefCov[rownames(coefCovR), colnames(coefCovR)] <- coefCovR
      
      for (i in 1:length(d)){
        if (d[[i]]$categorical) {
          coefCov[ rownames(d[[i]]$coefCov), 
                   colnames(d[[i]]$coefCov) ] <- d[[i]]$coefCov
        }
      }
    
    } # end restrictions
  
  } # end missing != listwise

  return(coefCov)
}