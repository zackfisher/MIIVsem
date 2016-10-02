#' @keywords internal
optimMIIV <- function(d, overid, data){
for (i in 1:length(d)){
  y  <- as.matrix( cbind(data[,d[[i]]$DVobs] ) )
  P  <- d[[i]]$IVobs
  k  <- overid + length(d[[i]]$IVobs) 
  
  if (k > length(d[[i]]$IV) ) {
    d[[i]]$NOTE <- paste("* Maximum number of MIIVs is less than requested degree of overidentification. See df for degree of overidentification.", sep="")
    d[[i]]$MSG <- "*"
  }
  
  else if (k <= length(d[[i]]$IV) ) {
    cd <- cor(data)
    cd <- cbind(cd[,P])
    cd[cd == 1] <- 0
    cd_m <- as.matrix(cd)
    
    cd_m <- cbind( apply(cd_m, 1, max) ) 
    ord_names <- rownames(cd_m)[order(cd_m, 
                                      decreasing=TRUE)][1:nrow(cd_m)]
    
    temp <- d[[i]]$IV
    temp <- temp[order(match(temp,ord_names))]
    d[[i]]$IV <- temp[1:k]
  }
}
return(d)
} 