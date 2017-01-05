#' @keywords internal
generateFormulas <- function(d, instruments){
  
  for (i in 1:length(d)){
    
      # Mikko: creating a formula based on a string is inefficient, and I do not see why we would
      # need to define the formulas here anyway because they would not be need when working
      # with covariance matrices.
      
      # Zack: I doubt it is inefficient, computationally.  Not the safest move though. 
      # The point about the utility of formulas for the covariance-based routines is a 
      # good one to keep in mind.  I'll separate the two functions for now. 
      
      d[[i]]$EqFormula <- stats::reformulate(d[[i]]$IVobs, d[[i]]$IVobs)
      d[[i]]$MIIVsFormula <- stats::reformulate(d[[i]]$IVobs, d[[i]]$IVobs)
   
  } 
  
  return(d)
  
}