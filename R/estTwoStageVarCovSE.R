#' estimate the SEs using two-stage ML
#'@keywords internal
estTwoStageVarCovSE <- function(g, v,eqns,pt){
  
      # fill the parTable 
      pt <- fillParTable(eqns, pt)
  
      # generate the model for variance covariance point 
      # estimates
      var.cov.model <- buildVarCovSyntax(pt)
      
      
      fit_vcov <- lavaan::lavaan(
        var.cov.model,
        sample.cov  = g$sample.cov,
        sample.mean = g$sample.mean,
        sample.nobs = g$sample.nobs,
        sample.cov.rescale = FALSE
      )
      

      pt_01 <- lavaan::parameterTable(fit_vcov)
      pt_01 <- pt_01[pt_01$op != "~*~" ,] 
      
      var.cov.params <- pt_01[pt_01$op == "~~" ,c("lhs","op","rhs","est")]
      var.cov.est <- var.cov.params$est
      names(var.cov.est) <- paste0(var.cov.params$lhs,var.cov.params$op,var.cov.params$rhs)
      
      #------------------------------------------------------------------------------#
      # Begin Theta2 Procedure
      #------------------------------------------------------------------------------#
      es <- lavaan::inspect(fit_vcov, what = "estimates")
      ef <- lavaan::inspect(fit_vcov, what = "free")
      LA <- as.matrix(es$lambda)
      TH <- as.matrix(es$theta)
      PS <- as.matrix(es$psi)
      BA<- as.matrix(es$beta) 
      
      # define vector
      LA.keep <- MIIVsem:::vec(LA)!=0 & MIIVsem:::vec(LA)!= 1
      LA.values <- MIIVsem:::vec(LA)[LA.keep]
      names(LA.values) <- MIIVsem:::vec(t(outer(colnames(LA), rownames(LA), fn3)))[LA.keep]
      
      TH.keep <- MIIVsem:::vec(ef$theta)!=0 
      TH.values <- MIIVsem:::vec(TH)[TH.keep]
      names(TH.values) <- MIIVsem:::vec(outer(colnames(ef$theta), rownames(ef$theta), fn2))[TH.keep]
      
      PS.keep <- MIIVsem:::vec(PS)!=0 
      PS.values <- MIIVsem:::vec(PS)[PS.keep]
      names(PS.values) <- MIIVsem:::vec(t(outer(colnames(PS), rownames(PS), fn2)))[PS.keep]
      
      BA.keep <- MIIVsem:::vec(BA)!=0 
      BA.values <- MIIVsem:::vec(BA)[BA.keep]
      names(BA.values) <- MIIVsem:::vec(outer(rownames(BA), colnames(PS), fn1))[BA.keep]
      
      
      # zf: commented out the line below bc we don't have access to mm$TH, which
      #     is essentially giving us the fixed and free params for the saturated
      #     model-impled covariance matrix. in the case of a correlation matrix 
      #     the diagonal would not be free. for now let's pretend we only have 
      #     continuous data
      lavsig <- lavaan::inspect(fit_vcov, what = "sigma")
      mat_ones <- matrix(0, nrow(lavsig), ncol(lavsig))
      mat_ones[lower.tri(mat_ones, diag = TRUE)] <- 1
      lavsig.free <- MIIVsem:::vec(mat_ones) != 0 
      
      EstSigmaTheta2 <-function(x, LA, TH, PS, BA,TH.k, PS.k, lavsig.free){
        TH[TH.k] <- x[1:sum(TH.k)]
        PS[PS.k] <- x[(sum(TH.k)+1):(sum(TH.k)+sum(PS.k))]
        SIG        <- LA %*% solve(diag(ncol(LA)) - BA) %*% PS %*% t(solve(diag(ncol(LA)) - BA)) %*% t(LA) + TH
        vechSigT   <- t(MIIVsem:::vec(SIG)[lavsig.free])
      }
      
      Delta2 <- numDeriv::jacobian(EstSigmaTheta2, c(TH.values, PS.values),LA=LA, TH=TH, PS=PS, BA=BA, TH.k=TH.keep, PS.k=PS.keep,lavsig.free=lavsig.free)
      names2 <- MIIVsem:::vec(outer(rownames(lavsig), colnames(lavsig), function(x, y) paste0(x, "~~", y)))[lavsig.free]
      rownames(Delta2) <- names2
      colnames(Delta2) <- names(c(TH.values, PS.values))
      
      EstSigmaTheta1 <-function(x, LA, TH, PS, BA, LA.k, BA.k,lavsig.free){
        LA[LA.k] <- x[1:sum(LA.k)]
        BA[BA.k] <- x[(sum(LA.k)+1):(sum(LA.k)+sum(BA.k))]
        SIG        <- LA %*% solve(diag(ncol(LA)) - BA) %*% PS %*% t(solve(diag(ncol(LA)) - BA)) %*% t(LA) + TH
        vechSigT   <- t(MIIVsem:::vec(SIG)[lavsig.free])
      }
      
      x <- c(LA.values, BA.values)
      Delta1 <- numDeriv::jacobian(EstSigmaTheta1, x,LA=LA, TH=TH, PS=PS, BA=BA,LA.k=LA.keep,BA.k=BA.keep,lavsig.free=lavsig.free)
      rownames(Delta1) <- names2
      colnames(Delta1) <- names(x)
      
      pos.db <- unlist(lapply(seq_along(rownames(Delta2)), function(ii){
        which(colnames(g$asymptotic.cov.sat) == rownames(Delta2)[ii] | colnames(g$asymptotic.cov.sat) == flip(rownames(Delta2)[ii]))
      }))
      wmat <- solve(g$asymptotic.cov.sat[pos.db,pos.db])
      H2 <- solve(t(Delta2) %*% wmat %*% Delta2)%*%t(Delta2) %*% wmat
      
      SigRho <- g$asymptotic.cov.sat
      SigRho <- SigRho[!grepl("\\|", colnames(SigRho)),!grepl("\\|", colnames(SigRho))]
      
      K <- matrix(0, nrow = length(names2) , ncol = length(c(LA.values,BA.values)) )
      rownames(K) <- names2
      colnames(K) <- unlist(lapply(eqns, function(x) { colnames(buildK(x, g$sample.cov))}))
      
      col_cnt <- 1
      
      for(jj in 1:length(eqns)){
        
        J <- buildK(eqns[[jj]], g$sample.cov, Svv_type = "arb")
        
        col_start <- col_cnt
        
        pos.db <- unlist(lapply(seq_along(rownames(J)), function(iii){
          which(rownames(K) == rownames(J)[iii] | rownames(K) == flip(rownames(J)[iii]))
        }))
        
        col_end   <- col_cnt + (ncol(J) - 1)
        K[pos.db, c(col_start:col_end)] <- J
        col_cnt <- col_end + 1
        
      }   
      
      K <- t(K)
      P1 <- H2 %*% (diag(nrow(Delta1)) - Delta1 %*% K)
      P2 <- t(diag(nrow(Delta1)) - Delta1 %*% K) %*% t(H2)
      
      
      pos.db <- unlist(lapply(seq_along(colnames(P1)), function(ii){
        which(rownames(SigRho) == colnames(P1)[ii] | rownames(SigRho) == flip(colnames(P1)[ii]))
      }))
      
      coefCov <- P1 %*% SigRho[pos.db, pos.db] %*% P2
      
      
      # # FIXME: remove suppress warnings.
      # # temporarily added suppress warnings here as 
      # # lavaan throws warnings re convergence 
      # # when iter.max is set to "none"
      # suppressWarnings( 
      #   var.cov.se.fit <- lavaan::lavaan(
      #     model = var.cov.se.model,
      #     sample.cov = g$sample.cov, 
      #     sample.mean = g$sample.mean,
      #     sample.nobs = g$sample.nobs,
      #     control=list(iter.max="none")
      #   )
      # )
      # 
      # DA <- unclass(lavaan::lavInspect(var.cov.se.fit, "delta"))
      # DA <- DA[c(grep("~~", rownames(DA)),grep("~1", rownames(DA))), ]
      # 
      # mi.cov   <- lavaan::inspect(var.cov.se.fit, "cov.ov")
      # mi.cov   <- mi.cov[rownames(g$sample.cov), colnames(g$sample.cov)]
      # mi.mean  <- lavaan::inspect(var.cov.se.fit, "mean.ov")
      # mi.mean  <- mi.mean[names(g$sample.mean)]
      # 
      # D        <- buildDuplication(length(g$sample.mean))
      # sigInv   <- solve(mi.cov)
      # meanDiff <- g$sample.mean - mi.mean
      # 
      # UL <- t(D) %*% 
      #   kronecker(
      #     sigInv, 
      #     sigInv %*% (g$sample.cov + tcrossprod(meanDiff)) %*% sigInv - .5*sigInv
      #   ) %*% D
      # 
      # LL <- kronecker(
      #   sigInv, 
      #   t(meanDiff) %*% sigInv 
      # ) %*% D
      # 
      # # Savalei (2010, Equation 5)
      # JB <- cbind(rbind(UL, LL), rbind(t(LL), sigInv))
      # 
      # ## Savalei (2010, Equation 7)
      # coefCov <- solve(t(DA) %*% JB %*% DA) %*% 
      #   t(DA) %*% JB %*% g$asymptotic.cov.sat %*% JB %*% 
      #   DA %*% solve(t(DA) %*% JB %*% DA)
      # 
      # param.names <- colnames(DA)
      # 
      # old.names <- param.names[grep("=~", param.names)]
      # 
      # if (length(old.names) > 0){
      #   
      #   new.names <- unlist(lapply(strsplit(old.names, "=~"), function (x){
      #     paste0(x[2],"~",x[1])
      #   }))
      #   
      #   for( i in 1:length(old.names)){
      #     param.names <- gsub(old.names[i], new.names[i], param.names)
      #   }
      #   
      # }
      # 
      # colnames(coefCov) <- rownames(coefCov) <- param.names
      
  return(coefCov)
}