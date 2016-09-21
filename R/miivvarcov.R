#' @keywords internal
miivvarcov <- function(model, estimator, dat, data, covariance){
  
  fit <- lavaan(model, 
                auto.fix.first = TRUE, 
                auto.var = TRUE, 
                auto.cov.lv.x = TRUE)
  
  lavsyntax <- lavExport(fit, export = FALSE)
  
  ls <- strsplit(lavsyntax, "\n") 
  
  for (i in 1:nrow(dat)){
    dv  <- dat[i,"DV"]
    iv  <- dat[i,"EV"]
    est <- dat[i,"Estimate"]
    
    for (m in 1:length(ls[[1]])){ 
      line <- ls[[1]][m] 
      if (grepl("=~", line)){
        line <- strsplit(line, "=~")
        
        if (iv == gsub(".*\\*","",line[[1]][1]) & 
            dv == gsub(".*\\*","",line[[1]][2])){
          ls[[1]][m] <- paste(iv, " =~ ", est, "*", dv,sep="")
        }
      }
      if (grepl(" ~ ", line)){
        line <- strsplit(line, "~")
        
        if (iv == gsub(".*\\*","",line[[1]][1]) & 
            dv == gsub(".*\\*","",line[[1]][2])){
          ls[[1]][m] <- paste(iv, " ~ ", est, "*", dv,sep="")
        }
      }
      if (grepl("==", line)){
        ls[[1]][m] <- NA
      }
    }
  }
  
  ls[[1]]  <- ls[[1]][!is.na(ls[[1]])]
  lavsyntax <- paste(unlist(ls), "\n", collapse="")
  lavsyntax <- gsub("NA", "1", lavsyntax)
  if (covariance == FALSE){
    fitcov   <- lavaan(lavsyntax, data = data, 
                       auto.fix.first = TRUE, estimator = estimator,
                       auto.var=TRUE, auto.cov.lv.x = TRUE)
  }
  if (covariance == TRUE){
    fitcov   <- lavaan(lavsyntax, sample.cov = cov, 
                       sample.nobs = N, auto.fix.first = TRUE, 
                       estimator = estimator, auto.var=TRUE, 
                       auto.cov.lv.x = TRUE)
  }
  pe  <-  parameterEstimates(fitcov)
  cpe <-  pe[pe$op == "~~" & pe$lhs != pe$rhs,]
  if (dim(cpe)[1] != 0){
    cpe$x <- paste(cpe$lhs, cpe$op, cpe$rhs, sep=" ")
    cpe <- cpe[,c("x","est", "se", "z", "pvalue")]
    colnames(cpe) <- c( "","Estimate", 
                        "StdErr", "Z", "P(|Z|)")
  }
  if (dim(cpe)[1] == 0){ cpe <- NULL }
  
  vpe <- pe[pe$op == "~~" & pe$lhs == pe$rhs,]
  vpe <- vpe[,c("lhs","est", "se", "z", "pvalue")]
  if (dim(vpe)[1] != 0){
    colnames(vpe) <- c( "","Estimate", "StdErr", "Z", "P(|Z|)")
  }
  if (dim(vpe)[1] == 0){ vpe <- NULL }
  vcov <- list(cov = cpe, var = vpe)
  
  return(list(vcov,lavsyntax))
}