#' estimate the 2SLS equation specific sigma squared
#'@keywords internal
est2SLSSigmaSq <- function(d, cov.mat){
  
    d <- lapply(d, function(eq) {
    
      if (eq$categorical){
      
        eq$sigma   <- NA
      
      } else {
      
        eq$sigma <- 
          (cov.mat[eq$DVobs, eq$DVobs] + 
              (t(eq$coefficients[-1]) %*%
          cov.mat[c(eq$IVobs), c(eq$IVobs)] %*% 
            eq$coefficients[-1]) -
            (2 * cov.mat[eq$DVobs, c(eq$IVobs)] %*% 
            eq$coefficients[-1]))
      
      }
      
      eq

    })

  return(d)
}