#' partial derivatives for PIV coef cov.
#'@keywords internal

dTheta.dSvv <- function(Svv, Svz, Svy, Svv_type = "cov", numeric = FALSE){
  
  U <- solve(t(Svz) %*% solve(Svv) %*% Svz)
  V <- t(Svz) %*% solve(Svv) %*% Svy
  
  # d vec [t(Svz) %*% solve(Svv) %*% Svz]
  # -------------------------------------
  # d vec  Svv 
  
  #d.Ua.Svv <- -1*kronecker(t(Svz)%*%t(solve(Svv)), t(Svz) %*%solve(Svv))
  
  # d vec solve[t(Svz) %*% solve(Svv) %*% Svz]
  # ------------------------------------------
  # d vec Svv      
  
  #d.U.Svv <- -1*kronecker(t(U), U) %*% d.Ua.Svv 
  
  
  # d vec [t(Svz) %*% solve(Svv) %*% Svy]
  # ------------------------------------------
  # d vec Svv    
  
  #d.V.Svv <- -1*kronecker(t(Svy)%*%t(solve(Svv)), t(Svz) %*%solve(Svv))
  
  # d vec [solve(U)%*%V]
  # ------------------------------------------
  # d vec Svv    
  
  #d.UV.Svv <- t(kronecker(V, diag(nrow(U)))) %*% d.U.Svv + 
  #            kronecker(diag(ncol(V)), U) %*% d.V.Svv
  
  
  # d vec [solve(U)%*%V]
  # ------------------------------------------
  # d vech Svv    
  
  #dTh.dSvv <- d.UV.Svv %*% makeD(ncol(Svv))
  
  # d vec [solve(U)%*%V]
  # ------------------------------------------
  # d vech Svv    
  
  k1 <- kronecker(t(V) %*% t(U) %*% t(Svz) %*% t(solve(Svv)), U %*% t(Svz) %*% solve(Svv))
  k2 <- kronecker(t(Svy) %*% t(solve(Svv)), U %*% t(Svz) %*% solve(Svv))
  
  # zf: commented out on 04/01/2019
  #
  if(Svv_type == "cov"){

    colnames.dTh.dSvv <- vech(outer(colnames(Svv), colnames(Svv), function(x, y) {  paste0(x, "~~", y) }))

    dTh.dSvv <- (k1-k2) %*% dupMat(nrow(Svv), ifelse(Svv_type=="cov", "s", "l"))
    
  } else if(Svv_type == "cor"){

    colnames.dTh.dSvv <- vecp(outer(colnames(Svv), colnames(Svv), function(x, y) {  paste0(x, "~~", y) }))

    dTh.dSvv <- (k1-k2) %*% dupMat(nrow(Svv), ifelse(Svv_type=="cov", "s", "l"))
    
  } else if (Svv_type == "arb"){
    
    # zf: added on 04/01/2019
    Svv.t <- Svv
    Svv.t[vec(Svv)==1]<-0
    K <- arbDupMat(Svv.t)
    dTh.dSvv <- (k1-k2) %*% K
    colnames.dTh.dSvv <- vec(outer(colnames(Svv), colnames(Svv), function(x, y) {  paste0(x, "~~", y) }))
    vecM <- c(MIIVsem:::vec(Svv.t))
    names.keep <- !(vecM==0 | duplicated(vecM))
    colnames.dTh.dSvv <- colnames.dTh.dSvv[names.keep]
    # end added on 04/01/2019
    
  }
  # end commented code
  

  colnames(dTh.dSvv) <- colnames.dTh.dSvv
  
  rownames.dTh.dSvv <- outer(colnames(Svy), colnames(Svz), function(x, y) {  paste0(x, "~", y) })
  
  rownames(dTh.dSvv) <- rownames.dTh.dSvv
  
  
  if(numeric){
    
   stop("numeric partial derivative not implemented.")
    
    dTh.dSvv.n <- function(vechSvv, Svz, Svy) {
      Svv <- revVech(vechSvv)
      solve(t(Svz) %*% solve(Svv) %*% Svz) %*% t(Svz) %*% solve(Svv) %*% Svy
    }
    
    dTh.dSvv.n <- numDeriv::jacobian(dTh.dSvv.n, vech(Svv),Svz = Svz,Svy = Svy)
    
    colnames(dTh.dSvv.n) <- colnames(dTh.dSvv)
    
    rownames(dTh.dSvv.n) <- rownames(dTh.dSvv)
    
    return(dTh.dSvv.n)
    
  } else {
    
    return(dTh.dSvv)
    
  }
  
  
}


##----------------------------------------------------------------------------
## Derivatives of theta w.r.t. Svz
#-----------------------------------------------------------------------------
dTheta.dSvz <- function(Svv, Svz, Svy, numeric = FALSE){
  
  U  <- solve(t(Svz) %*% solve(Svv) %*% Svz)
  V  <- t(Svz) %*% solve(Svv) %*% Svy
  TH <- U %*% V
  
  # d vec [solve(U)%*%V]
  # ------------------------------------------
  # d vec Svz    
  
  P3 <- kronecker(t(Svy) %*% t(solve(Svv)), U) %*% comMat(nrow(Svz), ncol(Svz)) 
  P1 <- kronecker(t(TH) %*% t(Svz) %*% t(solve(Svv)), U) %*% comMat(nrow(Svz), ncol(Svz))
  P2 <- kronecker(t(TH), U %*% t(Svz) %*% solve(Svv)) 
  dTh.dSvz <-  P3 - P1 - P2
  
  # P1 <- kronecker(t(Svy) %*% t(solve(Svv)), U) 
  # P3 <- kronecker(t(TH) %*% t(Svz) %*% t(solve(Svv)), U)
  # P2 <- kronecker(t(TH), U %*% t(Svz) %*% solve(Svv)) 
  # dTh.dSvz <-  (P3 - P1) %*%  comMat(nrow(Svz), ncol(Svz))  - P2
  
  # names
  
  colnames.dTh.dSvz <- vec(outer(rownames(Svz), colnames(Svz), function(x, y) {  paste0(x, "~~", y) }))
  
  colnames(dTh.dSvz) <- colnames.dTh.dSvz
  
  rownames.dTh.dSvz <- outer(colnames(Svy), colnames(Svz), function(x, y) {  paste0(x, "~", y) })
  
  rownames(dTh.dSvz) <- rownames.dTh.dSvz
  
  
  if(numeric){
    
    dTh.dSvz.n <- function(Svz, Svv, Svy) {
      solve(t(Svz) %*% solve(Svv) %*% Svz) %*% t(Svz) %*% solve(Svv) %*% Svy
    }
    
    dTh.dSvz.n <- numDeriv::jacobian(dTh.dSvz.n, Svz, Svv = Svv, Svy = Svy)
    
    colnames(dTh.dSvz.n) <- colnames.dTh.dSvz
    
    rownames(dTh.dSvz.n) <- rownames.dTh.dSvz
    
    return(dTh.dSvz.n)
    
  } else {
    
    return(dTh.dSvz)
    
  }
  
  
}


##----------------------------------------------------------------------------
## Derivatives of theta w.r.t. Svy
#-----------------------------------------------------------------------------
dTheta.dSvy <- function(Svv, Svz, Svy, numeric = FALSE){
  
  U <- solve(t(Svz) %*% solve(Svv) %*% Svz)
  
  # d vec solve[t(Svz) %*% solve(Svv) %*% Svz] %*% 
  #             t(Svz) %*% solve(Svv) %*% Svy
  # ----------------------------------------------
  # d vec  Svy
  
  dTh.dSvy <- U %*% t(Svz) %*% solve(Svv)
  
  # names
  
  colnames.dTh.dSvy <- vec(outer(rownames(Svy), colnames(Svy), function(x, y) {  paste0(x, "~~", y) }))
  
  colnames(dTh.dSvy) <- colnames.dTh.dSvy
  
  rownames.dTh.dSvy <- outer(colnames(Svy), colnames(Svz), function(x, y) {  paste0(x, "~", y) })
  
  rownames(dTh.dSvy) <- rownames.dTh.dSvy
  
  
  if(numeric){
    
    dTh.dSvy.n <- function(Svy, Svv, Svz) {
      solve(t(Svz) %*% solve(Svv) %*% Svz) %*% t(Svz) %*% solve(Svv) %*% Svy
    }
    
    dTh.dSvy.n <- numDeriv::jacobian(dTh.dSvy.n, Svy, Svv = Svv, Svz = Svz)
    
    colnames(dTh.dSvy.n) <- colnames.dTh.dSvy
    
    rownames(dTh.dSvy.n) <- rownames.dTh.dSvy
    
    return(dTh.dSvy.n)
    
  } else {
    
    return(dTh.dSvy)
    
  }
  
  
}

##----------------------------------------------------------------------------
## Derivatives of intercepts w.r.t. Svv
#-----------------------------------------------------------------------------
dIntercept.dSvv <- function(Svv, Svz, Svy, MuDV, MuIV, Svv_type = "cov", numeric = TRUE){
  
  if(Svv_type == "cov"){

    colnames.dInt.dSvv <- vech(outer(colnames(Svv), colnames(Svv), function(x, y) {  paste0(x, "~~", y) }))
    rownames.dInt.dSvv <- paste0(names(MuDV),"~1")
    
    dInt.dSvv.n <- function(vechSvv, Svz, Svy, MuDV, MuIV) {
      Svv <- revVech(vechSvv)
      B   <- solve(t(Svz) %*% solve(Svv) %*% Svz) %*% t(Svz) %*% solve(Svv) %*% Svy
      MuDV - MuIV%*%B
    }
    
    dInt.dSvv.n <- numDeriv::jacobian(dInt.dSvv.n, vech(Svv),Svz = Svz,Svy = Svy, MuDV=MuDV,MuIV=MuIV)
    
    colnames(dInt.dSvv.n) <- colnames.dInt.dSvv
    
    rownames(dInt.dSvv.n) <- rownames.dInt.dSvv
    
    return(dInt.dSvv.n)

  } else if(Svv_type == "cor"){

    colnames.dInt.dSvv <- vecp(outer(colnames(Svv), colnames(Svv), function(x, y) {  paste0(x, "~~", y) }))
    rownames.dInt.dSvv <- paste0(names(MuDV),"~1")

  } else if (Svv_type == "arb"){
    
      Svv.free <- Svv
      Svv.free[Svv.free == 1] <- 0
      k <- arbDupMat(Svv.free)
      
      x <- MIIVsem:::vech(Svv)
      x <- x[MIIVsem:::vech(Svv.free) != 0]
      set.to.one <- c(MIIVsem:::vec(Svv) == 1)
    
      dInt.dSvv.n <- function(x, Svz, Svy, MuDV, MuIV, k , set.to.one) {
        Svv <- MIIVsem:::revVec(k %*% t(t(x)))
        Svv[set.to.one]  <- 1
        B   <- solve(t(Svz) %*% solve(Svv) %*% Svz) %*% t(Svz) %*% solve(Svv) %*% Svy
        MuDV - MuIV%*%B
      }
    
      dInt.dSvv.n <- numDeriv::jacobian(dInt.dSvv.n, x,Svz = Svz,Svy = Svy, MuDV=MuDV,MuIV=MuIV, k=k , set.to.one=set.to.one)
    
      Svv.t <- Svv
      Svv.t[vec(Svv)==1]<-0
      cn <- vec(outer(colnames(Svv), colnames(Svv), function(x, y) {  paste0(x, "~~", y) }))
      vecM <- c(MIIVsem:::vec(Svv.t))
      names.keep <- !(vecM==0 | duplicated(vecM))
      colnames(dInt.dSvv.n) <- cn[names.keep]
      rownames(dInt.dSvv.n) <- paste0(names(MuDV),"~1")
      
      return(dInt.dSvv.n)

  }
  
}

##----------------------------------------------------------------------------
## Derivatives of intercept w.r.t. Svz
#-----------------------------------------------------------------------------
dIntercept.dSvz <- function(Svv, Svz, Svy, MuDV, MuIV, numeric = TRUE){
  
  colnames.dInt.dSvz <- vec(outer(rownames(Svz), colnames(Svz), function(x, y) {  paste0(x, "~~", y) }))
  rownames.dInt.dSvz <- paste0(names(MuDV),"~1")
  
  
  if(numeric){
    
    dInt.dSvz.n <- function(Svz, Svv, Svy, MuDV, MuIV) {
      B   <- solve(t(Svz) %*% solve(Svv) %*% Svz) %*% t(Svz) %*% solve(Svv) %*% Svy
      MuDV - MuIV%*%B
    }
    
    dInt.dSvz.n <- numDeriv::jacobian(dInt.dSvz.n, Svz, Svv = Svv, Svy = Svy, MuDV=MuDV,MuIV=MuIV)
    
    colnames(dInt.dSvz.n) <- colnames.dInt.dSvz
    
    rownames(dInt.dSvz.n) <- rownames.dInt.dSvz
    
    return(dInt.dSvz.n)
    
  } else {
    
    #return(dInt.dSvz)
    
  }
  
  
}



##----------------------------------------------------------------------------
## Derivatives of theta w.r.t. Svy
#-----------------------------------------------------------------------------
dIntercept.dSvy <- function(Svv, Svz, Svy, MuDV, MuIV, numeric = TRUE){
  
  colnames.dInt.dSvy <- vec(outer(rownames(Svy), colnames(Svy), function(x, y) {  paste0(x, "~~", y) }))
  rownames.dInt.dSvy <- paste0(names(MuDV),"~1")
  

  if(numeric){
    
    dInt.dSvy.n <- function(Svy, Svv, Svz, MuDV, MuIV) {
      B   <- solve(t(Svz) %*% solve(Svv) %*% Svz) %*% t(Svz) %*% solve(Svv) %*% Svy
      MuDV - MuIV%*%B
    }
    
    dInt.dSvy.n <- numDeriv::jacobian(dInt.dSvy.n, Svy, Svv = Svv, Svz = Svz, MuDV=MuDV,MuIV=MuIV)
    
    colnames(dInt.dSvy.n) <- colnames.dInt.dSvy
    
    rownames(dInt.dSvy.n) <- rownames.dInt.dSvy
    
    return(dInt.dSvy.n)
    
  } else {
    
    #return(dInt.dSvy)
    
  }
  
  
}