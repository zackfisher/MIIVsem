#' build model-implied correlation matrix
#' 
#' @param eqns 
#' @param pt 
#' 
#'@keywords internal
buildMICor <- function(eqns, pt){ 
  
  # Fill parTable with fixed regression coefficients.
  r    <- unlist(sapply(eqns,"[[", c("coefficients")))
  pt.z <- cbind(do.call(rbind, strsplit(names(r), "~")), r)
  pt.z <- pt.z[which(pt.z[,2]!="1"),] 

  for(i in 1:nrow(pt.z)){
    eq <- pt.z[i,]
    pt[pt$op == "=~" & pt$lhs %in% eq[2] & pt$rhs %in% eq[1], 
       c("free", "ustart")] <- c(0, as.numeric(eq[3]))
    pt[pt$op == "~"  & pt$lhs %in% eq[1] & pt$rhs %in% eq[2], 
       c("free", "ustart")] <- c(0, as.numeric(eq[3]))
  }
}

