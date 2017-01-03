#' @keywords internal
prepareRawData <- function(data) {
  
  callNoDots <- match.call(expand.dots = FALSE) 
  modelFrame <- callNoDots[c(1, match("data", names(callNoDots), 0))]
  modelFrame$na.action <- as.name("na.pass")
  modelFrame[[1]] <- as.name("model.frame")
  
  return( modelFrame )
}