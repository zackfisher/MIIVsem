#'@keywords internal
buildK <- function(eq, poly.mat,Svv_type = "cov"){
  
  Svv <-   poly.mat[eq[["MIIVs"]], eq[["MIIVs"]], drop = FALSE]
  Svy <-   poly.mat[eq[["MIIVs"]], eq[["DVobs"]], drop = FALSE]
  Svz <-   poly.mat[eq[["MIIVs"]], eq[["IVobs"]], drop = FALSE]
  
  dTheta.dSvv.d <- dTheta.dSvv(Svv, Svz, Svy, Svv_type)
  dTheta.dSvz.d <- dTheta.dSvz(Svv, Svz, Svy)
  dTheta.dSvy.d <- dTheta.dSvy(Svv, Svz, Svy)
  
  jac <- t(cbind(cbind(dTheta.dSvv.d, dTheta.dSvz.d), dTheta.dSvy.d))
  
  jac
  
}
