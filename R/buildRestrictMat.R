#' build the restriction matrix
#' 
#' Function returns a list \code{r} with information regarding
#' coefficient restrictions and estimation. The field \code{R} 
#' contains a restriction matrix for imposing cross-coefficient
#' restrictions. If no cross-coefficient restrictions are present
#' \code{R} is set to \code{NULL}. The vector \code{q} contains
#' a vector of lagrangean multipliers. The vector \code{constrained} 
#' gives the constrained coefficients in convenient form.
#'  
#' @param d A list containing equation information.
#' 
#' @keywords internal
buildRestrictMat <- function(d){
  
  # Reurns a vector of length equal to the 
  # number of coefficients estimated in the full
  # model. Coefficients without labels are listed
  # as NA. 
  coefIsLabeled <- unlist(lapply(d, function(eq){
    if (any(eq$categorical)) eq$Label else c(NA, eq$Label)
  }))
  
  # Labels may still be used without imposing restrictions.
  # Here we remove those labels. 
  coefIsRestricted <- sapply(coefIsLabeled, function(x) {
    if (!is.na(x) & (is.numeric(utils::type.convert(x, as.is=TRUE)) | 
        length(coefIsLabeled[coefIsLabeled %in% x]) > 1)){
      TRUE
    } else {
      FALSE
    } 
  })
  
  eq.restricted <- unlist(sapply(d, function (eq){
    ifelse(any(eq$Label %in% coefIsLabeled[coefIsRestricted]), TRUE, FALSE)
  }))
  names(eq.restricted) <- unlist(lapply(d, "[[", "DVobs"))
  
  
  # If there are no restricted coefficients...
  if (all(!coefIsRestricted)){
    
    R    <- NULL
    q    <- NULL
    cons <- NULL
  
  # ...or we have restricted coefficients
  } else {
    
    # Returns a vector containing the dependent
    # variable and independent variable pairs for
    # each coefficient in the model. If there
    # are any categorical variables in an equation
    # the intercept term is omitted from the equation.
    coefPairs <- unlist(lapply(d, function(eq){
      if (any(eq$categorical)) { 
        paste0(eq$DVobs,"_",eq$IVobs) 
      } else { 
        paste0(eq$DVobs,"_",c("(Intercept)",eq$IVobs)) 
      }
    }))
    


    # remove NAs from both vectors
    coef_i <- coefPairs[coefIsRestricted]
    labs_i <- coefIsLabeled[coefIsRestricted]
    
    cons <- c()
      
    for(l in unique(labs_i)){ 
        
      # if the label is numeric
      if (is.numeric(utils::type.convert(l, as.is=TRUE))){
          
        tmp  <- paste0(coef_i[labs_i == l],"=",labs_i[labs_i == l])
        cons <- c(cons, tmp)
          
      # the label is not numeric
      } else { 
          
        # choose the first contrained coefficient to 
        # be the pivot for all subsequent contraints 
        # in case there are more than 2 coefficients 
        # per label. If there is only one coefficient
        # per label, ignore it.
            
        primary   <- coef_i[labs_i == l][ 1L]
        secondary <- coef_i[labs_i == l][-1L]
            
        tmp  <- paste0(primary,"=",secondary)
        cons <- c(cons, tmp)

      }
    }
      
    R <- car::makeHypothesis(coefPairs, cons,rhs = NULL)
    if(is.null(dim(R))){ R  <- t( R) }
    q  <- R[,  ncol(R), drop = FALSE ]
    R  <- R[, -ncol(R), drop = FALSE ]
    rownames(R) <-  car::printHypothesis(R, q, colnames(R))
  }
  
  r <- list(R = R, q =q, constrained = cons, eq.restricted = eq.restricted)
  
  return(r)
}