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


derPvv <- function(Pvv, Pvz, Pvy){
  if(dim(Pvv)[1] == 1){
    deriv = NULL
    names = NULL
  } else {
    U1    <- solve(t(Pvz) %*% solve(Pvv) %*% Pvz)
    U2    <- t(Pvz) %*% solve(Pvv) %*% Pvy
    TERM1 <- kronecker(t(solve(Pvv)%*%Pvz%*%U1 %*%U2),U1 %*% t(Pvz)%*%solve(Pvv))
    TERM2 <- kronecker(t(solve(Pvv)%*%Pvy),U1 %*% t(Pvz)%*%solve(Pvv)) 
    deriv <- t((TERM1 - TERM2) %*% buildDuplicator(nrow(Pvv)))
    names <- t(utils::combn(colnames(Pvv), 2))
  }
  return(list(deriv = deriv, names = names))
}


derPvz <- function(Pvv, Pvz, Pvy){
  if(identical(Pvv,Pvz)){
    deriv = NULL
    names = NULL
  } else {
    U     <- solve(t(Pvz) %*% solve(Pvv) %*% Pvz)
    TH    <- U %*% t(Pvz) %*% solve(Pvv) %*% Pvy
    ER    <- (Pvy) - (Pvz %*% TH) 
    deriv <- t(kronecker(U, t(solve(Pvv) %*% ER)) - 
               kronecker(t(TH), U %*% t(Pvz) %*% solve(Pvv)))
    names <- as.matrix(expand.grid(rownames(Pvz),colnames(Pvz)))
  }
  return(list(deriv = deriv, names = names))
}


derPvy <- function(Pvv, Pvz, Pvy){
  deriv <- t(solve(t(Pvz) %*% solve(Pvv) %*% Pvz) %*% t(Pvz) %*% solve(Pvv))
  names <- as.matrix(expand.grid(rownames(Pvy),colnames(Pvy)))
  return(list(deriv = deriv, names = names))
}

buildCategoricalK <- function(eq, mat){
  
  eq$regDerivatives <- derRegPIV(
    mat[eq$MIIVs, eq$MIIVs, drop = FALSE], 
    mat[eq$MIIVs, eq$IVobs, drop = FALSE], 
    mat[eq$MIIVs, eq$DVobs, drop = FALSE]
  )
  
  # acmPos is a matrix containing the rho_{i,j} indices where i 
  # (column 1) and j (column 2) refer to the position of the 
  # variable (column 3) in the symmetric polychoric  
  # correlation matrix. 
  
  acmPos <- apply(t(utils::combn(colnames(mat), 2)), 2, function(x){
    match(x, colnames(mat))
  })
  acmPos <- cbind(acmPos, c(1:nrow(acmPos)))  
  
  reg.varID <- rbind(apply(eq$regDerivatives$names, 2, function(x){
    match(x, colnames(mat))
  }))
  
  # In reg.varID we need to order the rows so that col1 < col2
  Var1 <- NULL; Var2 <- NULL
  reg.varID <- transform(reg.varID, 
                         min = pmin(Var1, Var2), 
                         max = pmax(Var1, Var2))[,-c(1:2)]
  
  # merge changes the order so this needs to be corrected.
  reg.varID$order <- seq(1:nrow(reg.varID))
  posInfo <- merge(reg.varID, acmPos, by=c(1,2))
  posInfo <- posInfo[order(posInfo$order), ]
  
  rowK <- matrix(0, ncol(eq$regDerivatives$deriv), ncol(utils::combn(colnames(mat), 2)))
  rowK[,posInfo[,4]] <- t(eq$regDerivatives$deriv)
  rowK
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

buildDuplication <- function(x){
  mat <- diag(x)
  index <- seq(x*(x+1)/2)
  mat[ lower.tri( mat , TRUE ) ] <- index
  mat[ upper.tri( mat ) ] <- t( mat )[ upper.tri( mat ) ]
  outer(c(mat), index , function( x , y ) ifelse(x==y, 1, 0 ) )
}