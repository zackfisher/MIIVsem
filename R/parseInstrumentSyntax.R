#' Parse the instrument information
#' 
#' Check to determine if (1) the user-supplied dependent variables exist in the
#  set of estimating equations and (2) if the instruments provided by the user are valid MIIVs.
#' 
#' @param d A list containing equation information
#' @param instruments A list of user-supplied instruments. By defaults, \code{instruments} is set to \code{NULL} 
#' @param miiv.check Logical indicating whether the validity of MIIVs should be checked
#' 
#' @details 
#' \describe{
#'   \item{\code{instruments}}{ 
#'     Using the \code{instruments} option you can specify the MIIVs directly 
#'     for each equation in the model. To utilize this option you must first 
#'     define a list of instruments using the syntax displayed below. After 
#'     the list is defined, set the \code{instruments} argument equal to the 
#'     name of the list of MIIVs. Note, \code{instruments} are specified for 
#'     an equation, and not for a specific endogenous variable. Also, the
#'     \code{instruments} argument can be used to prune the estimated model.
#'     For example, if the \code{instruments} argument is provided only those
#'     equationss whose dependent variable is specified will be estimated.
#'   }
#'   \item{\code{miiv.check}}{ 
#'     The default is true.\code{miiv.check} can be used to turn off MIIV 
#'     validity cehcks, for example, if instruments are not in model.
#'   }
#' }
#'
#'@examples
#' 
#' \dontrun{
#' 
#'   #-----------------------------------------------------#
#'   # Examples:                                           #
#'   #  Ex.1 (Valid miivs provided)                        #
#'   #-----------------------------------------------------#
#' 
#'   model <- '
#'      Eta1 =~ y1 + y2  + y3  + y4  
#'      Eta2 =~ y5 + y6  + y7  + y8    
#'      Xi1  =~ x1 + x2 + x3 
#'
#'      Eta1 ~ Xi1  
#'      Eta2 ~ Xi1 
#'      Eta2 ~ Eta1 
#'
#'      y1   ~~ y5
#'      y2   ~~ y4
#'      y2   ~~ y6
#'      y3   ~~ y7
#'      y4   ~~ y8
#'      y6   ~~ y8 
#'   '
#'   
#'   good_instruments <- ' 
#'     y1 ~ x2 + x3                            
#'     y5 ~ y2 + y3 + y4 + x2                
#'     y2 ~ y3 + y7 + y8 + x2           
#'     y3 ~ y2 + y4 + y6 + y8        
#'     y4 ~ y3 + y6           
#'     y6 ~ y3 + y4 + y7 + x2            
#'     y7 ~ y2 + y4 + y6 + y8       
#'     y8 ~ y2 + y3 + y7 + x2          
#'     x2 ~ y1 + y5 + y2 + y3 + y4 + y6
#'     x3 ~ y1 + y5 + y2 + y3 + y4 + y6
#'   '
#'   
#'   bad_instruments <- ' 
#'     y1 ~ x2 + x3 + y4                         
#'     y5 ~ y2 + y3 + y4 + x2                
#'     y2 ~ y3 + y7 + y8 + x2           
#'     y3 ~ y2 + y4 + y6 + y8        
#'     y4 ~ y3 + y6           
#'     y6 ~ y3 + y4 + y7 + x2            
#'     y7 ~ y2 + y4 + y6 + y8       
#'     y8 ~ y2 + y3 + y7 + x2          
#'     x2 ~ y1 + y5 + y2 + y3 + y4 + y6
#'     x3 ~ y1 + y5 + y2 + y3 + y4 + y6
#'   '
#' 
#' }
#'   
#'@keywords internal
parseInstrumentSyntax <- function(d, instruments, miiv.check){
  
  if (!is.null(instruments) & miiv.check){
    
    mt   <- lavaan::lavParTable(instruments)
    mt   <- mt[mt$op == "~",]
    mts  <- split(mt[,c("lhs", "rhs")], mt$lhs)
  
    # Are all the dvs from the instrument list in the eqns list?
    if (length(setdiff(names(mts), lapply(d,"[[","DVobs"))) > 0) {
      
      stop("MIIVsem: dependent variables ", 
           paste0(setdiff(names(mts), lapply(d,"[[","DVobs")),collapse="," ),
           " not found in the set of estimating equations.")
      
    }
    
    # Are all the MIIVs from the instrument valid MIIVs?
    for (z in names(mts)){
      
      badMiivs <- setdiff(mts[[z]]$rhs,d[lapply(d,"[[","DVobs")==z][[1]]$MIIVs)
      
      if (length(badMiivs) > 0){
        
        stop("MIIVsem: The MIIV(s) ", paste0(badMiivs, collapse =", "),
             " specified for dependent variable ", z, 
             " are not valid instruments given the model syntax.")
        
      } else {
        d[lapply(d,"[[","DVobs")==z][[1]]$MIIVs <- mts[[z]]$rhs
      }
    }
    
    # Replace the pre
  }
  
  if (!is.null(instruments) & !miiv.check){
    
    mt   <- lavaan::lavParTable(instruments)
    mt   <- mt[mt$op == "~",]
    mts  <- split(mt[,c("lhs", "rhs")], mt$lhs)
    
    for (z in names(mts)){
      d[lapply(d,"[[","DVobs")==z][[1]]$MIIVs <- mts[[z]]$rhs
    }
  }
  
  return(d)
  
}