# parameters in test.stat: 
# 1. $sargan.statis classic(naive) Sargan.
# 2. $chi.sq, $mean.scale, $mean.var.adjust are adjust Sargan, mean-scaled statistic, mean-variance adjusted statistic, respectively.
# 3. $chi.sq.adj, $mean.scale.tilde, $mean.var.adjust.tilde are the omega tilde version of them, respectively.

over.stat.eq <- function(eq,g){
  # 1. Input
  # n <- nrow(data) # sample size
  # eq.num <- 1 # no. of example equation
  # iv <- c()
  # iv$eqns <- model1$eqn
  # iv.eq <- model1$eqn[[eq.num]] #information including DV, IV, MIIV
  n <- g$sample.nobs
  iv.eq <- eq
  
  ####      s_vz: s_vz used for estimation
  ####      s_vv: s_vv used for estimation
  ####      s_vu: s_vu used for estimation (u is y in the paper)
  # sample.covcor <- model1$sample.polychoric
  
  # also computed stat for continuous
  if (eq$categorical) {
      sample.covcor <- g$sample.polychoric
  } else {
    sample.covcor <- g$sample.cov
  }
  
  # fit <- lavaan::lavCor( object = reisenzein1986, output = "fit", missing = "listwise", estimator = "two.step", se = "standard", ov.names.x = NULL, ordered = paste("Z", c(5 : 7, 11 : 13), sep = "") )
  # sample.covcor <- unclass(lavaan::inspect(fit, "cov.ov"))
  s_vz <- sample.covcor[ iv.eq$MIIVs, iv.eq$IVobs, drop = FALSE ]
  s_vv <- sample.covcor[ iv.eq$MIIVs, iv.eq$MIIVs, drop = FALSE ]
  s_vu <- sample.covcor[ iv.eq$MIIVs, iv.eq$DVobs, drop = FALSE ]
  
  inv_s_vzzzzv <- solve( t(s_vz) %*% solve(s_vv) %*% s_vz )
  ####      theta: estimated theta_j^(1)
  # matrix(eq[["coefficients"]][2])
  
  # # 2025/10 ----
  iv.coef <-  eq[["coefficients"]]
  iv.coef <- iv.coef[!names(iv.coef) %in% paste0(eq$DVlat,"~1")]  # omit constant value
  theta <-  matrix(iv.coef) # omit constant value
  colnames(theta) <- eq$DVobs
  rownames(theta) <- eq$IVobs
  
  # following is correct only when continuous variable
  # theta <-  matrix(eq[["coefficients"]][-1]) # omit constant value
  # colnames(theta) <- eq$DVobs
  # rownames(theta) <- eq$IVobs
  
  # following is correct only when no parameters were fixed
  # theta <- inv_s_vzzzzv %*% t(s_vz) %*% solve(s_vv) %*% s_vu ## We estimate theta here to ensure the correct set of IVs is used.
  # 2025/10 ----
  
  ####      sigma.names: rownames/colnames of the upsilon matrix (asymptotic covariance matrix)
  # sigma.names <- model1$asymptotic.cov %>% colnames
  sigma.names <- colnames(g$asymptotic.cov)
  
  ####      ordered: names of ordinal variables
  # ordered <- model1$ordered
  
  ####      upsilon: asymptotic covariance matrix upsilon
  upsilon <- n*g$asymptotic.cov
  
  # 2. Calculate omega
  # computes d[g_j(sigma)] / d[sigma]
  dgsigma_dsigma <- matrix( 0, nrow(s_vu), length(sigma.names) )
  colnames(dgsigma_dsigma) <- sigma.names
  
  # i1 <- 1
  # i2 <- 1
  ####  Derivatives with respect to the elements in s_vu
  row_name_svu <- rownames(s_vu)
  col_name_svu <- colnames(s_vu)
  for( i1 in 1 : nrow(s_vu) ){
    for( i2 in 1 : ncol(s_vu) ){
      col.pos <- which( sigma.names == paste(row_name_svu[i1], '~~', col_name_svu[i2], sep = '') )
      if( length(col.pos) == 0 ){
        col.pos <- which( sigma.names == paste(col_name_svu[i2], '~~', row_name_svu[i1], sep = '') )
      }
      ####  Initiate the derivative of s_vu with respect to sigma
      der_s_vu <- matrix( 0, nrow(s_vu), ncol(s_vu) )
      der_s_vu[i1, i2] <- 1
      ####  Get the gradient
      dgsigma_dsigma[, col.pos] <- der_s_vu
      
    }
  }
  
  ####  Derivatives with respect to the elements in s_vz
  row_name_svz <- rownames(s_vz)
  col_name_svz <- colnames(s_vz)
  for( i1 in 1 : nrow(s_vz) ){
    for( i2 in 1 : ncol(s_vz) ){
      col.pos <- which( sigma.names == paste(row_name_svz[i1], '~~', col_name_svz[i2], sep = '') )
      if( length(col.pos) == 0 ){
        col.pos <- which( sigma.names == paste(col_name_svz[i2], '~~', row_name_svz[i1], sep = '') )
      }
      ####  Initiate the derivative of s_vu with respect to sigma
      der_s_vz <- matrix( 0, nrow(s_vz), ncol(s_vz) )
      der_s_vz[i1, i2] <- 1
      ####  Get the gradient
      dgsigma_dsigma[, col.pos] <- -1.0 * der_s_vz %*% theta
      
    }
  }
  
  dgsigma_dsigma
  
  # computes omega
  omega <- dgsigma_dsigma %*% upsilon %*% t(dgsigma_dsigma)
  
  ## 3.1 Adjusted Sargan chi square statistic
  no.iv <- nrow(s_vu)
  no.endog <- ncol(s_vz)
  
  ####  Start to compute test statistics
  test.stat <- list() ## It is initiated as a list, but will be unlist() in the end.
  resid.miiv <- s_vu - s_vz %*% theta
  
  ####  Omega matrix and its square root matrix
  eigen.omega <- eigen(omega)
  omega.sqrt <- eigen.omega$vectors %*% diag(sqrt(eigen.omega$values)) %*% t(eigen.omega$vectors)
  inv.omega.sqrt <- solve(omega.sqrt)
  
  ####  Qmat and the MP inverse of Q * t(Q)
  ####  We used the ginv() function. It can also be computed using eigen-decomposition
  qmat <- diag(1,no.iv) - inv.omega.sqrt %*% s_vz %*% inv_s_vzzzzv %*% t(s_vz) %*% solve(s_vv) %*% omega.sqrt
  gmat <- MASS::ginv( qmat %*% t(qmat) )
  
  ####  The chi-square distributed test statistic, Eq.(19). It is also the generalized Wald statistic
  test.stat$chi.sq <- c( n * t(resid.miiv) %*% inv.omega.sqrt %*% gmat %*% inv.omega.sqrt %*% resid.miiv )
  test.stat$chi.sq.df <- no.iv - no.endog
  test.stat$chi.sq
  test.stat$chi.sq.df
  
  ## 3.2 Omega tilda  and corresponding Adjusted Sargan chi square statistic
  omega.adj <- omega - ( s_vu %*% t(s_vu) - s_vu %*% t(theta) %*% t(s_vz) - s_vz %*% theta %*% t(s_vu) + s_vz %*% theta %*% t(theta) %*% t(s_vz) )
  eigen.omega.adj <- eigen(omega.adj)
  if( all(eigen.omega.adj$values >= 0) ){
  omega.sqrt.adj <- eigen.omega.adj$vectors %*% diag(sqrt(eigen.omega.adj$values)) %*% t(eigen.omega.adj$vectors)
  inv.omega.sqrt.adj <- solve(omega.sqrt.adj)
  qmat.adj <- diag(1,no.iv) - inv.omega.sqrt.adj %*% s_vz %*% inv_s_vzzzzv %*% t(s_vz) %*% solve(s_vv) %*% omega.sqrt.adj
  gmat.adj <- MASS::ginv( qmat.adj %*% t(qmat.adj) )
  test.stat$chi.sq.adj <- c( n * t(resid.miiv) %*% inv.omega.sqrt.adj %*% gmat.adj %*% inv.omega.sqrt.adj %*% resid.miiv )
  test.stat$chi.sq.adj
  } else {
    test.stat$chi.sq.adj <- NA
  }
  
  ## 3.3  (Naive) Sargan test 
  s_uu <- sample.covcor[ iv.eq$DVobs, iv.eq$DVobs, drop=FALSE ]
  s_zu <- sample.covcor[ iv.eq$IVobs, iv.eq$DVobs, drop=FALSE ]
  s_zz <- sample.covcor[ iv.eq$IVobs, iv.eq$IVobs, drop=FALSE ]
  
  phi2_j <- c( s_uu - t(theta) %*% s_zu - t(s_zu) %*% theta + t(theta) %*% s_zz %*% theta )
  sargan.stat <- n * c( t(resid.miiv) %*% solve(s_vv) %*% resid.miiv ) / phi2_j
  test.stat$sargan.stat <- sargan.stat
  test.stat$sargan.stat
  ## 3.4 Satorra-Bentler statistics 
  
  pi.mat <- t(qmat) %*% omega.sqrt %*% solve(s_vv) %*% omega.sqrt %*% qmat / phi2_j
  # equivalent to pi.mat <- omega.sqrt %*% ( solve(s_vv) - solve(s_vv) %*% s_vz %*% inv_s_vzzzzv %*% t(s_vz) %*% solve(s_vv) ) %*% omega.sqrt / phi2_j
  test.stat$mean.scale <- (no.iv - no.endog) * sargan.stat / sum(diag(pi.mat))
  test.stat$mean.var.adjust <- sum(diag(pi.mat)) * sargan.stat / sum(diag(pi.mat %*% pi.mat))
  test.stat$mean.var.df <- (sum(diag(pi.mat)) ^ 2) / sum(diag(pi.mat %*% pi.mat))
  test.stat$mean.scale
  test.stat$mean.var.adjust
  test.stat$mean.var.df
  
  ## 3.5  Satorra-Bentler statistics using tilde omega
  if( all(eigen.omega.adj$values >= 0) ){
  pi.mat.adj <- t(qmat.adj) %*% omega.sqrt.adj %*% solve(s_vv) %*% omega.sqrt.adj %*% qmat.adj / phi2_j
  test.stat$mean.scale.tilde <- (no.iv - no.endog) * sargan.stat / sum(diag(pi.mat.adj))
  test.stat$mean.var.adjust.tilde <- sum(diag(pi.mat.adj)) * sargan.stat/sum(diag(pi.mat.adj %*% pi.mat.adj))
  test.stat$mean.var.df.tilde <- (sum(diag(pi.mat.adj)) ^ 2) / sum(diag(pi.mat.adj %*% pi.mat.adj))
  test.stat$mean.scale.tilde
  test.stat$mean.var.adjust.tilde
  test.stat$mean.var.df.tilde
  } else {
    test.stat$mean.scale.tilde <- NA
    test.stat$mean.var.adjust.tilde <- NA
    test.stat$mean.var.df.tilde <- NA
  }
  
  return(test.stat)
  
}