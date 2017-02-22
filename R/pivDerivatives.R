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
  numDerPvv     <- function(vechPvv, Pvz, Pvy){
    matdim <- ceiling(sqrt(2*length(vechPvv)))
    Pvv    <- matrix(0, nrow = matdim, ncol = matdim)
    Pvv[col(Pvv) < row(Pvv)] <- vechPvv
    Pvv <- Pvv + t(Pvv)
    diag(Pvv) <- 1
    solve(t(Pvz)%*%solve(Pvv)%*%Pvz) %*% t(Pvz)%*%solve(Pvv) %*% Pvy
  }
  deriv <- t(numDeriv::jacobian(numDerPvv, Pvv[lower.tri(Pvv)], Pvz = Pvz, Pvy = Pvy))
  names <- t(combn(colnames(Pvv), 2))
  return(list(deriv = deriv, names = names))
}

# Partial derivative of theta with respect to Pvz
derPvz <- function(Pvv, Pvz, Pvy){
  numDerPvz     <- function(Pvz, Pvv, Pvy){
    U  <- solve(t(Pvz) %*% solve(Pvv) %*% Pvz)
    U %*% (t(Pvz)%*% solve(Pvv) %*% Pvy)
  }
  deriv <- t(numDeriv::jacobian(numDerPvz, Pvz, Pvv = Pvv, Pvy = Pvy))
  names <- as.matrix(expand.grid(rownames(Pvz),colnames(Pvz)))
  return(list(deriv = deriv, names = names))
}

# Partial derivative of theta with respect to Pvy
derPvy <- function(Pvv, Pvz, Pvy){
  numDerPvy     <- function(Pvy, Pvz, Pvv){
    solve(t(Pvz)%*%solve(Pvv)%*%Pvz) %*% t(Pvz)%*%solve(Pvv) %*% Pvy
  }
  deriv <- t(numDeriv::jacobian(numDerPvy, Pvy, Pvz = Pvz, Pvv = Pvv))
  names <- as.matrix(expand.grid(rownames(Pvy),colnames(Pvy)))
  return(list(deriv = deriv, names = names))
}


buildKmatrix <- function(d, pcr){
  eq <- d[[1]]
  # acmPos is a matrix containing the rho_{i,j} indices where i 
  # (column 1) and j (column 2) refer to the position of the 
  # variable (column 3) in the symmetric polychoric  
  # correlation matrix. 
  acmPos <- apply(t(combn(colnames(pcr), 2)), 2, function(x){
    match(x, colnames(pcr))
  })
  acmPos <- cbind(acmPos, c(1:nrow(acmPos)))  
  # View(cbind(acmPos, colnames(acov)))
  krows <- lapply(d, function(eq){
    
    reg.varID <- apply(eq$regDerivatives$names, 2, function(x){
      match(x, colnames(pcr))
    })
    
    # In reg.varID we need to order the rows so that col1 < col2
    reg.varID <- transform(reg.varID, 
                           min = pmin(Var1, Var2), 
                           max = pmax(Var1, Var2))[,-c(1:2)]
    
    # merge changes the order so this needs to be corrected.
    reg.varID$order <- seq(1:nrow(reg.varID))
    posInfo <- merge(reg.varID, acmPos, by=c(1,2))
    posInfo <- posInfo[order(posInfo$order), ]
    
    rowK <- rep(0, 1/2 * nrow(pcr)* (nrow(pcr)-1))
    rowK[posInfo[,4]] <- eq$regDerivatives$deriv[,1]
    rowK
  })
  return(do.call("rbind", krows))
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