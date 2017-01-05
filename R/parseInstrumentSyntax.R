#' @keywords internal
parseInstrumentSyntax <- function(d, instruments){
  
  
  if (is.null(instruments)){
    for (i in 1:length(d)){
      d[[i]]$MIIVsUsed <- d[[i]]$MIIVs
    }
  } else if (class(instruments) == "character"){
    
    mt   <- lavParTable(instruments)
    mt   <- mt[mt$op == "~",]
    mts  <- split(mt[,c("lhs", "rhs")], mt$lhs)
  
    # Are all the dvs from the instrument list in the eqns list?
    if (length(setdiff(names(mts), lapply(d,"[[","DVobs"))) > 0){
      
      stop("MIIVsem: The dependent variable ", 
           setdiff(names(mts), lapply(d,"[[","DVobs")),
           " was not found in the set of estimating equations.")
      
    }
    
    # Are all the MIIVs from the instrument valid MIIVs?
    for (z in names(mts)){
      
      badMiivs <- setdiff(mts[[z]]$rhs,d[lapply(d,"[[","DVobs")==z][[1]]$MIIVs)
      
      if (length(badMiivs) > 0){
        
        stop("MIIVsem: The MIIV(s) ", badMiivs,
             " specified for dependent variable ", z, 
             " are not valid instruments given the model syntax.")
        
      }
    }
  }
  
  return(d)
  
}