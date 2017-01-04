#' @keywords internal
generateFormulas <- function(d, instruments){
  
  
  if (is.null(instruments)){
    for (i in 1:length(d)){
      d[[i]]$MIIVsUsed <- d[[i]]$MIIVs
      
      # Mikko: creating a formula based on a string is inefficient, and I do not see why we would
      # need to define the formulas here anyway because they would not be need when working
      # with covariance matrices.
      
      d[[i]]$EqFormula <- stats::as.formula(c(paste0(d[[i]]$DVobs, "~"), 
                                        paste0(d[[i]]$IVobs,collapse = "+")))
      d[[i]]$MIIVsFormula <- stats::as.formula(c(paste0("~"), 
                                           paste0(d[[i]]$MIIVs, collapse = "+")))
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
    
    for (i in 1:length(d)){
      
      d[[i]]$MIIVsUsed  <- d[[i]]$MIIVs
      d[[i]]$EqFormula  <- stats::as.formula(c(paste0(d[[i]]$DVobs, "~"), 
                                        paste0(d[[i]]$IVobs,collapse = "+")))
      
      if (d[[i]]$DVobs %in% names(mts)){
        
        tmpDV <- d[[i]]$DVobs
        d[lapply(d,"[[","DVobs")==tmpDV][[1]]$MIIVsUsed <- mts[[tmpDV]]$rhs
        
      }
      
      d[[i]]$MIIVsFormula  <- stats::as.formula(c(paste0("~"), 
                                           paste0(d[[i]]$MIIVsUsed, collapse = "+")))
      
    }
    
    
  }
  
  return(d)
  
}