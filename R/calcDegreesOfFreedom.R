#' @keywords internal
calcDegreesOfFreedom <- function(d, restrictions, numObsEq_i, numCoef_i){
  
  
  numObsAll    <- sum(unlist(numObsEq_i)) 
  numCoefUnr   <- sum(unlist(numCoef_i))
  numCoefRes_i <- numCoef_i
  
  if(!is.null(restrictions)) {
    R            <- restrictions$R
    numCoefRes   <- numCoefUnr - nrow(R)
    
    for(j in 1:nrow(R)) {
      for(i in 1:length(d)) {  
        # Identify restrictions that are NOT cross-equation
        colIx1 <- (1+sum(numCoef_i[1:i])-numCoef_i[i])
        colIx2 <- (sum(numCoef_i[1:i]))
        if( sum(R[j, colIx1:colIx2]^2) == sum(R[j,]^2)) {
          numCoefRes_i[i] <- numCoefRes_i[i] - 1
        }
      }
    }
  } else {
    numCoefRes_i <- numCoef_i
    numCoefRes   <- numCoefUnr 
  }
  
  df <- unlist(numObsEq_i) - numCoefRes_i   
  return(list(df = df,  numCoefRes= numCoefRes,  numCoefRes_i= numCoefRes_i))
}