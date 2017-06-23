#' Select MIIVs based on different criteria
#'  
#' @keywords internal
optimMIIV <- function(d, overid, sample.cov, sample.mean, sample.nobs){
  
  enum.miivs <- function(x, k) {
    
    if( k > length(x) ) {
      
      stop('overid > number of available MIIVs.')
      
    }
    
    if(choose(length(x), k)==1){
      
      list(as.vector(combn(x, k)))
    
    } else {
      
      cbn <- combn(x, k)
      lapply(seq(ncol(cbn)), function(i) cbn[,i])
      
    }
  }
  
  d <- lapply(d, function(eq){
    
    IVobs     <- eq$IVobs
    
    miiv.list <- enum.miivs(eq$MIIVs, length(eq$IVobs)+overid)
    
    stat.list <- unlist(lapply(miiv.list, function(MIIVs) {
      
      mean(sheasRSq(IVobs, MIIVs, sample.cov, sample.mean, sample.nobs)[,1])
      
    }))
    
    eq$MIIVs <- miiv.list[[which.max(stat.list)]]
    
    eq
    
  })
  
  
return(d)
  
} 