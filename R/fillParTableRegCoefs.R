#'@keywords internal
function <- estVarCovar(eqns, pt){
  # Fill parTable with fixed regression coefficients.
  r <- do.call(c,sapply(d,"[[", c("coefficients")))
  z <- cbind(do.call(rbind, strsplit(names(r), "~")), r)
  z <- z[which(z[,2]!="1"),] # cuts down on nrow

  for(i in 1:nrow(z)){
    eq <- z[i,]
    pt$ustart[pt$op == "=~" & pt$lhs %in% eq[2] & pt$rhs %in% eq[1]] <- as.numeric(eq[3])
    pt$ustart[pt$op == "~" & pt$lhs %in% eq[1] & pt$rhs %in% eq[2]] <- as.numeric(eq[3])
  }
  
  return(pt)
}

