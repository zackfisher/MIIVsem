#'@keywords internal

# Partial derivative of theta with respect to regression params
derRegPIV <- function(Pvv, Pvz, Pvy){

  dPvv   <- derPvv(Pvv, Pvz, Pvy)
  dPvz   <- derPvz(Pvv, Pvz, Pvy)
  dPvy   <- derPvy(Pvv, Pvz, Pvy)
  
  return(list(
    deriv = rbind(dPvv$deriv, dPvz$deriv, dPvy$deriv),
    names = rbind(dPvv$names, dPvz$names, dPvy$names))
  )
}

# Partial derivative of theta with respect to Pvv
derPvv <- function(Pvv, Pvz, Pvy){
  U1    <- solve(t(Pvz) %*% solve(Pvv) %*% Pvz)
  U2    <- t(Pvz) %*% solve(Pvv) %*% Pvy
  TERM1 <- kronecker(t(solve(Pvv)%*%Pvz%*%U1 %*%U2),U1 %*% t(Pvz)%*%solve(Pvv))
  TERM2 <- kronecker(t(solve(Pvv)%*%Pvy),U1 %*% t(Pvz)%*%solve(Pvv)) 
  
  deriv <- t((TERM1 - TERM2) %*% buildDuplicator(nrow(Pvv)))
  names <- t(combn(colnames(Pvv), 2))
  
  return(list(deriv = deriv, names = names))
}

# Partial derivative of theta with respect to Pvz
derPvz <- function(Pvv, Pvz, Pvy){
  U     <- solve(t(Pvz) %*% solve(Pvv) %*% Pvz)
  TH    <- U %*% t(Pvz) %*% solve(Pvv) %*% Pvy
  ER    <- (Pvy) - (Pvz %*% TH) 
  deriv <- t(kronecker(U, t(solve(Pvv) %*% ER)) - 
               kronecker(t(TH), U %*% t(Pvz) %*% solve(Pvv)))
  names <- as.matrix(expand.grid(rownames(Pvz),colnames(Pvz)))
  return(list(deriv = deriv, names = names))
}

# Partial derivative of theta with respect to Pvy
derPvy <- function(Pvv, Pvz, Pvy){
  deriv <- t(solve(t(Pvz) %*% solve(Pvv) %*% Pvz) %*% t(Pvz) %*% solve(Pvv))
  names <- as.matrix(expand.grid(rownames(Pvy),colnames(Pvy)))
  return(list(deriv = deriv, names = names))
}


buildKmatrix <- function(d, pcr){
  
  # replace variable names with their position in the dataframe 
  eq <- d[[1]]
  reg.varID <- apply(eq$regDerivatives$names, 2, function(x){
    match(x, colnames(pcr))
  })
  
}

vecp <- function( x ){return( t( t( x[lower.tri(x)] ) ) )}

vec  <- function( x ){ return( t( t( as.vector( x ) ) ) )}

buildDuplicator <- function(p){
  D <- matrix(0, nrow = p^2, ncol = .5*(p*(p-1)))
  I <- diag(rep(1, p))
  for (i in 1:p) {
    for (j in 1:p) {  
      if (i > j){
        
        vij <- rep(0, times = .5*(p*(p-1)))
        vij[(j-1)*p + i - .5*(j*(j+1))] <- 1
        
        eij <-  I[i, ] %o% I[j, ]
        eji <-  I[j, ] %o% I[i, ]
        tij <-  vec(eij + eji)
        
        D <- D + (tij %*% vij)
      }
    }
  }
  D
}