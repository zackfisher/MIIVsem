#' Prune MIIVs based on some specified criterion.
#'  
#' @keywords internal
pruneExcessMIIVs <- function(d, 
                             overid.degree, 
                             overid.method, 
                             data = NULL, 
                             sample.cov = NULL, 
                             sample.mean = NULL, 
                             sample.nobs = NULL){
  
  if (!overid.method %in% c("partial.R2","partial.r2","minimum.eigen", "minimum.eigenvalue", "random")){
          stop(paste0("miive: overid.method ", overid.method, " not supported."))
  }
  
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
    
    if(overid.method == "random"){
      
      eq$MIIVs <- sample(eq$MIIVs, (length(eq$IVobs)+ overid.degree), replace=FALSE)
      
    } else {
    
      miiv.list <- enum.miivs(eq$MIIVs, length(eq$IVobs)+overid.degree)
      
      stat.list <- unlist(lapply(miiv.list, function(MIIVs) {
        
        
        if (overid.method == "partial.R2" | overid.method == "partial.r2"){
          mean(sheasRSq(IVobs, MIIVs, sample.cov, sample.mean, sample.nobs)[,1])
        } else if (overid.method == "minimum.eigenvalue" | overid.method == "minimum.eigen"){
          minEigenStat(IVobs, MIIVs, data)
        }
        
      }))
      
      eq$MIIVs <- miiv.list[[which.max(stat.list)]]
      
    }
    
    eq
    
  })
  
return(d)
  
} 