#' @keywords internal 
processMissing <- function(Y_i, Z_i, V_i){
  
  noNA <- matrix(NA, nrow = length(Y_i[[1]]), ncol = length(Y_i))
  
  # An element in the noNA matrix is TRUE if...
  #  for the ith observation in the jth equation
  #    (1) the ith element in the jth DV is not missing 
  #    (2) there are no missing values in the Z_j matrix for observation i
  #    (3) there are no missing values in the X_j matrix for observation i
  for(i in 1:length(Y_i)){ # i <- 1
    noNA[,i] <- !is.na(Y_i[[i]]) & 
      rowSums(is.na(Z_i[[i]])) == 0 & 
      rowSums(is.na(V_i[[i]])) == 0
  }
  
  # Is the system unbalanced 
  noNA_j <- rowSums(!noNA) == 0
  
  for(i in 1:ncol(noNA)) {
    if(any(noNA[!noNA_j, i])) {
      stop( "MIIVsem: The system of equations is not balanced." )
    }
  }
  return(noNA)
}