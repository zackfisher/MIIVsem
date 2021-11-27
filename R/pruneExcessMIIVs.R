#' Prune MIIVs based on some specified criterion.
#'  
#' @keywords internal
pruneExcessMIIVs <- function(d, 
                             overid.degree, 
                             overid.method, 
                             data = NULL, 
                             sample.cov = NULL, 
                             sample.polychoric = NULL,
                             sample.mean = NULL, 
                             sample.nobs = NULL,
                             overid.boot = FALSE,
                             cat.vars    = NULL){
  
  
  if (!overid.method %in% c("random",
                            "stepwise.R2")){
    
          stop(paste0("miive: overid.method ", overid.method, " not supported."))
  }
  

  # if we are not in a bootstrap replication
  # save the original MIIV set as OrigMIIVs
  if(!overid.boot) {
    d <- lapply(d, function(eq){
      eq$OrigMIIVs <- eq$MIIVs
      eq
    })
  # or we are in a bootstrap loop and want to
  # reselect from the original MIIV set
  } else {
    d <- lapply(d, function(eq){
      eq$MIIVs <- eq$OrigMIIVs
      eq
    })
  }
  
  enum.miivs <- function(x, k) {
    
    if( k > length(x) ) {
      
      stop('overid > number of available MIIVs.')
      
    }
    
    if(choose(length(x), k)==1){
      
      list(as.vector(utils::combn(x, k)))
    
    } else {
      
      cbn <- utils::combn(x, k)
      lapply(seq(ncol(cbn)), function(i) cbn[,i])
      
    }
  }

  if (overid.method == "stepwise.R2"){
    
   d <- lapply(d, function(eq){
      
     exo.included <- intersect(eq$IVobs, eq$MIIVs)
     end.included <- setdiff(eq$IVobs, eq$MIIVs)
     
     eqn.all.miivs <- setdiff(eq$MIIVs, exo.included)
     eqn.miivs     <- c()
     eqn.R2        <- c()
     eqn.num.miivs <- length(end.included) + overid.degree
     
     if (eqn.num.miivs >= length(eq$MIIVs) | length(eqn.all.miivs) < eqn.num.miivs){
       
       eqn.miivs <- eq$MIIVs
       
     } else {
     
       if(is.null(sample.polychoric)){
         
         sample.cor <- stats::cov2cor(sample.cov)
         
       } else {
         
         #### ZF 2018-12-28
         #### is this always a covariance matrix? No.
         sample.cor <- stats::cov2cor(sample.polychoric)
         
       }
        
       #--------------------------------------------------# 
       # Obtain for each explanatory variable (EV) the MIIV 
       # with the largest squared pearson correlation 
       # coefficient  (R^2). Once a MIIV is selected for a
       # given EV it is no longer eligible to be selected 
       # for subsequent EVs.
       #--------------------------------------------------# 
       for (i in 1:length(end.included)){
         cor.row       <- sample.cor[end.included[i], eqn.all.miivs, drop = F]^2
         eqn.R2        <- c(eqn.R2, cor.row[which.max(cor.row)])
         eqn.miivs     <- c(eqn.miivs, colnames(cor.row)[which.max(cor.row)])
         names(eqn.R2) <- c(names(eqn.R2)[-length(eqn.R2)], end.included[i])
         eqn.all.miivs <- setdiff(eqn.all.miivs, eqn.miivs)
       }
       
       #--------------------------------------------------# 
       # While we still need for MIIVs...
       #   Select a new MIIV which results in the largest
       #   stepwise increase in R^2 for the EV with the 
       #   smallest R^2.
       #--------------------------------------------------# 
       while (length(eqn.miivs) < eqn.num.miivs){
         
         IVobs <- names(which.min(eqn.R2))
         
         row.R2 <- rbind(lapply(eqn.all.miivs, function(x){
           b    <- solve(sample.cor[c(eqn.miivs,x),c(eqn.miivs,x)], sample.cor[c(eqn.miivs,x), IVobs ])
           cors <- sample.cor[IVobs,c(eqn.miivs,x)]
           rsq  <- sum(b*cors)
           rsq
         }))
         
         # replaced on 12.27.2018 to remove dependence on raw data
         #
         # row.R2 <- rbind(lapply(eqn.all.miivs, function(x){
         #   formula <- stats::reformulate(c(eqn.miivs,x), response=IVobs)
         #   rsq <- summary(lm(formula, data = as.data.frame(data)))$r.squared
         #   rsq
         # }))
         
         colnames(row.R2) <- eqn.all.miivs
         eqn.R2[names(eqn.R2) %in% IVobs]    <- as.numeric(row.R2[which.max(row.R2)])
         eqn.miivs        <- c(eqn.miivs, colnames(row.R2)[which.max(row.R2)])
         eqn.all.miivs    <- setdiff(eqn.all.miivs, eqn.miivs)
  
       }
     }
     
     eq$MIIVs <- eqn.miivs
     eq
      
    })

  } else {
  
    d <- lapply(d, function(eq){
      
      IVobs     <- eq$IVobs
      
      if(overid.method == "random"){
        
        if(length(eq$MIIVs) > (length(eq$IVobs)+ overid.degree)){
          eq$MIIVs <- sample(eq$MIIVs, (length(eq$IVobs)+ overid.degree), replace=FALSE)
        } 
        
      } else if (overid.method == "lasso"){
        
        if(length(eq$MIIVs) > (length(eq$IVobs)+ overid.degree)){
          # exo.included <- intersect(eq$IVobs, eq$MIIVs)
          # end.included <- setdiff(eq$IVobs, eq$MIIVs)
          # eqn.all.miivs <- setdiff(eq$MIIVs, exo.included)
          # eqn.num.miivs <- length(end.included) + overid.degree
          # eq$MIIVs <- adaptiveLasso(eq$DVobs, end.included, eqn.all.miivs, data, eqn.num.miivs) 
        }
      } else {
      
        miiv.list <- enum.miivs(eq$MIIVs, length(eq$IVobs)+overid.degree)
        
        stat.list <- unlist(lapply(miiv.list, function(MIIVs) {
          
          
          # if (overid.method == "partial.R2" | overid.method == "partial.r2"){
          #   mean(sheasRSq(IVobs, MIIVs, sample.cov, sample.mean, sample.nobs)[,1])
          # } else if (overid.method == "minimum.eigenvalue" | overid.method == "minimum.eigen"){
          #   minEigenStat(IVobs, MIIVs, data)
          # }
          
        }))
        
        eq$MIIVs <- miiv.list[[which.max(stat.list)]]
        
      }
      
      eq
      
    })
  }
  return(d)
  
} 