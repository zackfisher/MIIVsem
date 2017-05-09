library("MIIVsem")

 #obj <- readRDS("~/Documents/pkgs/MIIVsem/MIIVsem/tests/testthat/rds/ex02_polr_coefcov.rds")

  model <- '
    Eta1 =~ y1 + l1*y2  + l2*y3  + l3*y4
    Eta2 =~ y5 + l1*y6  + l2*y7  + l3*y8
    Xi1  =~ x1 + x2 + .5*x3
    Eta1 ~ l4*Xi1
    Eta2 ~ l4*Xi1
    Eta2 ~ l5*Eta1

    y1   ~~ y5
    y2   ~~ y4
    y2   ~~ y6
    y3   ~~ y7
    y4   ~~ y8
    y6   ~~ y8
  '

  fit <- miive(model, bollen1989a)
  
  
 #-------------------------------------------------------# 
  context("ex02: poldemo (r) coefficients correct")
 #-------------------------------------------------------# 
  
  expect_equal_to_reference(
    fit$coefficients[sort(names(fit$coefficients))], 
    "rds/ex02_polr_coef.rds"
  )
  
 #-------------------------------------------------------# 
  context("ex02: poldemo (r) coefficient covariance matrix correct")
 #-------------------------------------------------------# 
  dim_names <- sort(names(fit$coefficients))
  
  expect_equal_to_reference(
    fit$coefCov[dim_names, dim_names], 
    "rds/ex02_polr_coefcov.rds"
  )

 #-------------------------------------------------------#  
  #context("ex02: poldemo (r) residual covariance matrix correct")
 #-------------------------------------------------------# 
  
  #dim_names <- sort(colnames(fit$residCov))
  
  #expect_equal_to_reference(
  #  fit$residCov[dim_names, dim_names], 
  #  "rds/ex02_polr_residcov.rds"
  #)
  
  #-------------------------------------------------------#  
   context("ex02: poldemo (r) sargan test statistics correct")
  #-------------------------------------------------------# 
  
  sargan <- unlist(lapply(fit$eqn, "[", c("sargan")))
  names(sargan) <- unlist(lapply(fit$eqn, "[", c("DVlat")))
  sargan   <- sargan[sort(names(sargan))]
  
  expect_equal_to_reference(
    sargan, 
    "rds/ex02_polr_sargan.rds"
  )
  
  #-------------------------------------------------------#  
  context("ex02: poldemo (r) sargan df correct")
  #-------------------------------------------------------# 
  
  sargandf <- unlist(lapply(fit$eqn, "[", c("sargan.df")))
  names(sargandf) <- unlist(lapply(fit$eqn, "[", c("DVlat")))
  sargandf   <- sargandf[sort(names(sargandf))]
  
  expect_equal_to_reference(
    sargandf, 
    "rds/ex02_polr_sargandf.rds"
  )
  
  #-------------------------------------------------------#  
  context("ex02: poldemo (r) sargan p-values correct")
  #-------------------------------------------------------# 
  
  sarganp <- unlist(lapply(fit$eqn, "[", c("sargan.p")))
  names(sarganp) <- unlist(lapply(fit$eqn, "[", c("DVlat")))
  sarganp   <- sarganp[sort(names(sarganp))]
  
  expect_equal_to_reference(
    sarganp, 
    "rds/ex02_polr_sarganp.rds"
  )

  