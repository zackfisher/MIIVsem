#' Calculate the residual covariance matrix
#' 
#'@keywords internal
estResidualCov <- function(d, sample.cov){
  
    # Calculate the residual covariance matrix for the full system of
    # equations based on covariance matrix input.
    
    dvs <- unlist(lapply(d, "[[", "DVobs"))
    B   <- diag(length(dvs))
    colnames(B) <- rownames(B) <- dvs
    idx <- do.call("rbind",lapply(d, function(eq){
      cbind(eq$DVobs, eq$IVobs, eq$coefficients[-1])
    } ))
    idx <- idx[idx[,2] %in% dvs, ,drop = FALSE]
    B[idx[,2:1, drop = FALSE]] <- -1*as.numeric(idx[,3])
    
    residCov <- t(B) %*% sample.cov[dvs,dvs] %*% B
     B %*% sample.cov[dvs,dvs] %*% t(B)
    
    colnames(residCov) <- rownames(residCov) <- unlist(lapply(d,"[[","DVlat"))
    residCov
    
  
}
