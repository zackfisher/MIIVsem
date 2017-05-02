#' Estimate the variance and covariance parameters
#' 
#' @param data
#' @param d
#' @param d.un
#' @param pt
#' @param ordered
#' 
#'@keywords internal
estVarCovar <- function(data, d, d.un, pt, ordered){
  
  v <- list()

  if(length(d.un) > 0){
    stop(paste("MIIVsem: variance covariance esimtation not",
                "allowed in the presence of underidentified",
                "equations."))
  }
    
  v$coefficients <- unlist(apply(
    lavaan::parameterEstimates(
      lavaan::sem(
        createModelSyntax(d, pt),
        data,
        ordered = ordered
      )), 1, function(x){
        if (x["op"] == "~~"){
          coefi <- x["est"]
          names(coefi) <- paste0(x["lhs"],x["op"],x["rhs"])
          coefi
        }
      }
    )
  )
  
  storage.mode(v$coefficients) <- "numeric"
  
  
  
  return(v)
}