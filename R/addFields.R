#' @keywords internal
addFields <- function(d, factorNames){
  
  #
  # Add logical "categorical" to list d indicating whether or not
  # the equation contains a categorical variable.
  # 
  if ( length(factorNames) > 0 ) { # add categorical label begin
    d <- lapply(d, function(eq){
      if(length(intersect(factorNames, c(eq$DVobs, eq$IVobs, eq$MIIVs)) > 0)){
        eq$categorical <- TRUE
      } else {
        eq$categorical <- FALSE
      }
      eq
    })
  } else {
    d <- lapply(d, function(eq){ eq$categorical <- FALSE; eq })
  } # add categorical label end
  
  
  
  labs_i <- do.call("rbind",(lapply(d, function(eq){
    cbind(eq$DVobs, eq$IVobs, eq$Label)
  })))
  
  # restricted dependent variables
  resDVs <- labs_i[duplicated(labs_i[,3,drop=FALSE]) & !is.na(labs_i[,3,drop=FALSE]),1]
  
  # add restricted field to equations list
  d <- lapply(d, function(eq){ 
    if (eq$DVobs %in% resDVs)  eq$restricted <- TRUE else eq$restricted <- FALSE
    eq 
  })
  
  return(d)
  
}