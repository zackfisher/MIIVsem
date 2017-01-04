#' @keywords internal
generateFormulas <- function(d, instruments){
  
  
  if (is.null(instruments)){
    for (i in 1:length(d)){
      d[[i]]$MIIVsUsed <- d[[i]]$MIIVs
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
  
  
#     # In the eqns list (d) add a 'MIIVsUsed' field containing the 
#     # the valid MIIVs to be used in estimation. Initially, this is
#     # all valid MIIVs. 
#     for (i in 1:length(d)){
#     
#       d[[i]]$MIIVsUsed <- d[[i]]$MIIVs
#       
#     }
# 
#     # Generate a list of formulas for the estimating equations.
#     asEqForm <- function(x){
#       stats::as.formula(c(paste0(x$DVobs, "~"), paste0(x$IVobs, collapse = "+")))
#     }
#     
#     eqFormulas <- lapply(d, asEqForm)
#     
# 
#   if (is.null(instruments)) {
#     
#     miivFormulas <- lapply(d, function(x) stats::as.formula(c(paste0("~",x$MIIVs))) )
#     
#   } else if (class(instruments) == "character"){
#     
#     mt   <- lavParTable(instruments)
#     mt   <- mt[mt$op == "~",]
#     
#     mts <- split(mt[,c("lhs", "rhs")], mt$lhs)
#     
#     asMiivForm <- function(x){
#       stats::as.formula(c(paste0("~"), paste0(x$rhs, collapse = "+")))
#     }
#     
#     miivFormulas <- lapply(mts, asMiivForm)
#   
#     # Perform some basic checks:
#     dvsFromInstruments <- names(miivFormulas)
#     dvsFromMiivsSearch <- unlist(lapply(d, "[", c("DVobs")), use.names = FALSE)
#     
#     # Are all the dvs from the instrument list in the eqns list?
#     if (length(setdiff(dvsFromInstruments, dvsFromMiivsSearch)) > 0){
#       
#       stop("MIIVsem: The dependent variable ", 
#            setdiff(dvsFromInstruments, dvsFromMiivsSearch),
#            " was not found in the set of estimating equations.")
#       
#     }
#     
# 
#     # Are all the MIIVs from the instrument valid MIIVs?
#     for (z in dvsFromInstruments){
#       
#       badMiivs <- setdiff(mts[[z]]$rhs,d[lapply(d,"[[","DVobs")==z][[1]]$MIIVs)
#       
#       if (length(badMiivs) > 0){
#         
#         stop("MIIVsem: The MIIV(s) ", badMiivs,
#              " specified for dependent variable ", z, 
#              " are not valid instruments given the model syntax.")
#         
#       }
#     }
#     
#     # Check if there are enough user-specified MIIVs to identify 
#     # the equations.
#     #
#     # Should we stop if there are more explanatory variables in an
#     # equation than instruments? For now, yes. 
# 
#     for (z in dvsFromInstruments){
#       
#       if (length(d[lapply(d,"[[","DVobs")==z]$IVobs) > length(mts[[z]]$rhs)){
#         
#         stop("MIIVsem: The number of explanatory variables ", 
#              " in the ", z, " equation exceed the ",
#              " number of user specified instruments, ",
#              length(mts[[z]]$rhs),   ".")
#         
#       }
#     }
#     
#     ## Edit the eqns list (d) for the specified equations and add
#     ## a 'MIIVsUsed' field in 'd' that contains the MIIVs to be 
#     ## used in the estimating equations.
#     
#     for (i in 1:length(d)){
#       
#       if (d[[i]]$DVobs %in% names(mts)){
#         
#         tmpDV <- names(mts)[which(d[[i]]$DVobs %in% names(mts))]
#         d[lapply(d,"[[","DVobs")==tmpDV]$MIIVsUsed <- mts[[tmpDV]]$rhs
#         
#       }
#       
#     }
#   }
# }