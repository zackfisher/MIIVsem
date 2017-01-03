#' @keywords internal
# returns NULL if there are no restrictions,
# otherwise returns a list containing the 'R' matrix
# and 'q' vector, as well as a vector 'cons' of 
# the constrained coefficients.
buildRestrictMat <- function(d){
  
  # vector containing dv_iv coefficient label
  coef_i <- unlist(lapply(d,  function(x) paste0(x$DVobs,"_",c("(Intercept)",x$IVobs))))
  
  
  # vector containing labels (add an NA term to correspond with intercepts)
  labs_i <- unlist(lapply(lapply(d, "[[", "Label"), function(x) c(NA, x)))
  
  # remove NAs from both vectors
  coef_j <- coef_i[!is.na(labs_i)]
  labs_i <- labs_i[!is.na(labs_i)]
  
  cons <- c()
  
  if(length(labs_i) > 0){
    for(l in unique(labs_i)){ # l <- "l1"
      # if the label is numeric
      if (is.numeric(type.convert(l, as.is=TRUE))){
        tmp  <- paste0(coef_j[labs_i == l],"=",labs_i[labs_i == l])
        cons <- c(cons, tmp)
      } else { # the label is not numeric
        if (length(labs_i[labs_i == l]) > 1L){
          # choose the first contrained coefficient to 
          # be the pivot for all subsequent contraints 
          # in case there are more than 2 coefficients 
          # per label.
          primary   <- coef_j[labs_i == l][ 1L]
          secondary <- coef_j[labs_i == l][-1L]
          tmp  <- paste0(primary,"=",secondary)
          cons <- c(cons, tmp)
        }
      }
    }
    
    R <- car::makeHypothesis(coef_i, cons,rhs = NULL)
    if(is.null(dim(R))){ R  <- t( R) }
    q  <- R[,  ncol(R), drop = FALSE ]
    R  <- R[, -ncol(R), drop = FALSE ]
    rownames(R) <-  car::printHypothesis(R, q, colnames(R))
    restrictList <- list(R = R, q = q, cons = cons)
    
  } else {
    
    # restrictList <- list(R = NULL, q = NULL, cons = NULL)
    restrictList <- NULL
    
  }
  return(restrictList)
}