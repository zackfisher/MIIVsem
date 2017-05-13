#' estimate the SEs using two-stage ML
#'@keywords internal
estTwoStageML <- function(g,v,eqns,pt){
  
  #-------------------------------------------------------#
  # two stage ML
  #-------------------------------------------------------#
      pt <- fillParTable(eqns, pt, v)
      
      var.cov.se.model <- buildMissingSyntax(pt)
      
      # FIXME: remove suppress warnings.
      # temporarily added suppress warnings here as 
      # lavaan throws warnings re convergence 
      # when iter.max is set to "none"
      suppressWarnings( 
        var.cov.se.fit <- lavaan::lavaan(
          model = var.cov.se.model,
          sample.cov = g$sample.cov, 
          sample.mean = g$sample.mean,
          sample.nobs = g$sample.nobs,
          control=list(iter.max="none")
        )
      )
      
      DA <- unclass(lavaan::lavInspect(var.cov.se.fit, "delta"))
      DA <- DA[c(grep("~~", rownames(DA)),grep("~1", rownames(DA))), ]
      
      mi.cov   <- lavaan::inspect(var.cov.se.fit, "cov.ov")
      mi.cov   <- mi.cov[rownames(g$sample.cov), colnames(g$sample.cov)]
      mi.mean  <- lavaan::inspect(var.cov.se.fit, "mean.ov")
      mi.mean  <- mi.mean[names(g$sample.mean)]
      
      D        <- buildDuplication(length(g$sample.mean))
      sigInv   <- solve(mi.cov)
      meanDiff <- g$sample.mean - mi.mean
      
      UL <- t(D) %*% 
        kronecker(
          sigInv, 
          sigInv %*% (g$sample.cov + tcrossprod(meanDiff)) %*% sigInv - .5*sigInv
        ) %*% D
      
      LL <- kronecker(
        sigInv, 
        t(meanDiff) %*% sigInv 
      ) %*% D
      
      JB <- cbind(rbind(UL, LL), rbind(t(LL), sigInv))

      ## Savalei (2010, Equation 7)
      coefCov <- solve(t(DA) %*% JB %*% DA) %*% 
        t(DA) %*% JB %*% g$asymptotic.cov.sat %*% JB %*% 
        DA %*% solve(t(DA) %*% JB %*% DA)
      
      param.names <- colnames(DA)
      
      old.names <- param.names[grep("=~", param.names)]
      
      if (length(old.names) > 0){
        
        new.names <- unlist(lapply(strsplit(old.names, "=~"), function (x){
          paste0(x[2],"~",x[1])
        }))
        
        for( i in 1:length(old.names)){
          param.names <- gsub(old.names[i], new.names[i], param.names)
        }
        
      }

      colnames(coefCov) <- rownames(coefCov) <- param.names
      
  return(coefCov)
}