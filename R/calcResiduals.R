#' @keywords internal
calcResiduals <- function(b, Y_i, Z_i, numCoef_i){
  
  b_i <- split(
    b, rep(1:length(numCoef_i), numCoef_i)
  )
  
  resids_i <- list()
  
  for(i in 1:length(Y_i)){
    
    resids_i[[i]] <- Y_i[[i]] - Z_i[[i]] %*% b_i[[i]]
    
  }
  
  return(resids_i)
  
}