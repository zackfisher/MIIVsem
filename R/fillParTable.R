#' fill the parameter table
#'@keywords internal
fillParTable <- function(eqns, pt, v = NULL){
  
    r <- unlist(lapply(eqns,"[[", c("coefficients")))
    z <- cbind(do.call(rbind, strsplit(names(r), "~")), r)
    
    for(i in 1:nrow(z)){
      
      eq <- z[i,]
      
      pt[
        pt$op == "=~" & pt$lhs %in% eq[2] & pt$rhs %in% eq[1], 
        c("ustart")
        ] <- c(as.numeric(eq[3]))
      
      pt[
        pt$op == "~" & pt$lhs %in% eq[1] & pt$rhs %in% eq[2], 
        c("ustart")
        ] <- c(as.numeric(eq[3]))
      
      # estimated intercepts
      pt[
        pt$op == "~1" & pt$lhs %in% eq[1] & eq[2] == "1", 
        c("ustart")
        ] <- c(as.numeric(eq[3]))
      
    }
    
    tmpMarkers <- pt[pt$op == "=~",]$rhs[which(!duplicated(pt[pt$op == "=~",]$lhs))]
    
    # now fix scaling indicator intercepts to zero. 
    pt[pt$op == "~1" & pt$lhs %in% tmpMarkers,  c("free","ustart")] <- c(0,0)
    
    # lavvan fixes exogenous LV means to zero, we need to free them. 
    latVars <- unique(pt$lhs[pt$op=="=~"])
    latEndVars <- unique(c(
      pt$rhs[pt$op == "=~" & pt$rhs %in% latVars], 
      pt$lhs[pt$op == "~"  & pt$lhs %in% latVars]
    )) 
    latExoVars <- setdiff(latVars, latEndVars)
    
    pt[pt$op == "~1" & pt$lhs %in% latExoVars,  c("ustart")] <- NA
    
    # free any non-zero latent variable regression intercepts
    pt[pt$op == "~1" & !is.na(pt$ustart) & as.numeric(pt$ustart)!= 0 & 
         pt$lhs %in%  pt[pt$op == "=~",]$lhs[
           which(!duplicated(pt[pt$op == "=~",]$lhs))
           ],  c("free")
       ] <- c(1)
    
    # reset labels for any latent variable regression coefficients
    pt[pt$op == "~1" & !is.na(pt$ustart) & as.numeric(pt$ustart)!= 0 & 
         pt$lhs %in%  pt[pt$op == "=~",]$lhs[
           which(!duplicated(pt[pt$op == "=~",]$lhs))
           ],  c("mlabel")
       ] <- NA
  
  if(!is.null(v)){
    
    z.v <- cbind(
      do.call("rbind", strsplit(names(v$coefficients), "~")), 
      v$coefficients
    )[,-2]
    
    
    for(i in 1:nrow(z.v)){
      eq <- z.v[i,]
      pt[pt$op == "~~" & pt$lhs %in% eq[1] & pt$rhs %in% eq[2], 
         c("ustart")
         ] <- c(as.numeric(eq[3]))
    }
    
  }
  return(pt)
  
}