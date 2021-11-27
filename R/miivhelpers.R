# Miscellanneous helper functions
# A helper function for binding a list of matrices by their names
#@param x a list of matrices
#@param empty value to be used for empty cells
#'@keywords internal
miivhelper.bindmatrices <- function(x, empty){
  
  coln <- unique(unlist(lapply(x,colnames)))
  rown <- unique(unlist(lapply(x,rownames)))
  
  mat <- matrix(empty,length(rown), length(coln),
                dimnames = list(rown,coln))
  
  for(i in seq_along(x)){
    mat[rownames(x[[i]]), colnames(x[[i]])] <- x[[i]]
  }
  
  mat
}


buildBlockDiag <- function(d, mat, row, col, inv = FALSE, pcr){
 
  res <- lavaan::lav_matrix_bdiag(
    if (inv){
      lapply(d, function(eq){
        solve(mat[c("1",eq[[row]]), c("1",eq[[col]]), drop = FALSE]) 
      })
    } else if (pcr) {
      lapply(d, function(eq){
        mat[eq[[row]], eq[[col]], drop = FALSE] 
      })
    } else {
      lapply(d, function(eq){
        mat[c("1",eq[[row]]), c("1",eq[[col]]), drop = FALSE] 
      })
    })
}

createModelSyntax <- function(eqns, pt){ # eqns <- results$eqn
  # Fill parTable with fixed regression coefficients.
  r <- unlist(lapply(eqns,"[[", c("coefficients")))
  z <- cbind(do.call(rbind, strsplit(names(r), "~")), r)
  #z <- z[which(z[,2]!="1"),] # cuts down on nrow
  
  for(i in 1:nrow(z)){
    eq <- z[i,]
    pt[pt$op == "=~" & pt$lhs %in% eq[2] & pt$rhs %in% eq[1], c("free", "ustart")] <- c(0, as.numeric(eq[3]))
    pt[pt$op == "~" & pt$lhs %in% eq[1] & pt$rhs %in% eq[2], c("free", "ustart")] <- c(0, as.numeric(eq[3]))
  }
  
  
  mod <- paste(apply(pt, 1, function(eq){
    ifelse(
      !is.na(eq["ustart"]),
      paste0(eq["lhs"], eq["op"], eq["ustart"], "*",eq["rhs"], "\n"),
      paste0(eq["lhs"], eq["op"],eq["rhs"], "\n")
    )
  }), collapse = "")
  
  mod.int <- paste(apply(z, 1, function(eq){
    if (eq[2] == "1"){
      paste0(eq[1], "~", eq[3], "*1", "\n")
    } else {
      paste0("")
    }
  }), collapse = "")

  mod <- paste0(mod, mod.int)
  
  return(mod)
}


makeName <- function(names, prefix = NULL) {
  W <- 14
  if(is.null(prefix)) {
    prefix <- rep("", length(names))
  }
  multiB <- FALSE
  if(any(nchar(names) != nchar(names, "bytes"))){
    multiB <- TRUE
  }
  if(!multiB) {
    names <- abbreviate(names, minlength = W, strict = TRUE)
  } else {
    names <- sprintf(paste("%-", W, "s", sep=""), names)
  }
  char.format <- paste("%3s%-", W, "s", sep = "")
  sprintf(char.format, prefix, names)
}


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


getIndicatorR2s <- function(indicators, data){
  return(unlist(lapply(seq_along(indicators), function(i){
    model <- stats::as.formula(paste(indicators[i],"~",paste(setdiff(indicators,indicators[i]), collapse="+") ))
    r2   <- summary(stats::lm(model, as.data.frame(data)))$r.squared
    names(r2) <- indicators[i]
    r2
  })))
}



# vec-type operators ====
vec   <- function(x) t(t(as.vector(x))) 
vech  <- function(x) t(t(x[!upper.tri(x)])) 
vecp  <- function(x) t(t(x[!upper.tri(x, diag = TRUE)]))
vecd  <- function(x) t(t(diag(x)))

# (reverse) vec-type operators ====

revVec <- function(x, ncol, nrow){
  
  if (length(x)==1)  { return(x) }
  
  d <- sqrt(length(x))
  
  if (missing(ncol) | missing(nrow)){
    ncol <- nrow <- d
    if (round(d) != d){
      stop("MIIVsem: Dimensions needed for non-square matrices.")
    }
  }
  
  revvecx <- matrix(0, nrow = nrow, ncol = ncol)
  
  for (j in 1:ncol){
    revvecx[,j] <- x[c(1:nrow) + (j-1)*nrow]
  }
  return(revvecx)
}

revVech <- function(x){
  
  if (length(x)==1) { return(x) }
  
  d <- (-1 + sqrt(8*length(x) + 1))/2
  
  if (round(d) != d)
    stop("MIIVsem: Matrix is not square.")
  
  revvechx <- matrix(0, nrow=d, ncol=d)
  
  for (j in 1:d){
    revvechx[j:d,j] <- x[1:(d-j+1)+ (j-1)*(d - 1/2*(j-2))]
  }
  
  revvechx <- revvechx + t(revvechx) - diag(diag(revvechx))
  
  return(revvechx)
}

revVecp <- function(x){
  
  d <- (1 + sqrt(8*length(x) + 1))/2
  
  if (round(d) != d){
    stop("MIIVsem: Matrix is not square.")
  }
  revvecp <- matrix(0, nrow=d, ncol=d)
  revvecp[lower.tri(revvecp , diag=FALSE)] <- x
  diag(revvecp) <- 1
  revvecp <- revvecp + t(revvecp) - diag(diag(revvecp))
  
  return(revvecp)
}

# basis functions ====

# DEF'N from Magnus & Neudecker (1980, p. 423): The unit
#   vector e_i, i,...,n, is the ith column of the identity
#   matrix I_n, or a vector length n with one in its 
#   ith position and zeros elsewhere.

make_e_basis <- function(i, n) {
  
  replace(numeric(n), i, 1)
  
}

# DEF'N from Magnus & Neudecker (1980, p. 423): The unit
#   vector of length [1/2 * n * (n + 1)] where there is
#   a one in the [(j-1) * n + i - .5*j*(j+1)] position
#   and zeros elsewhere 1 <= j <= i <= n.

make_u_basis <- function(i,j,n) {
  
  replace(numeric(.5*n*(n+1)), (j-1)*n + i - .5*j*(j-1), 1)
  
}

# DEF'N from Neudecker (1983, p. 274): The unit
#   vector of length [1/2 * n * (n - 1)] where there is
#   a one in the [(j-1) * n + i - .5*j*(j+1)] position
#   and zeros elsewhere 1 <= j < i <= n.

make_v_basis <- function(i,j,n) {
  
  replace(numeric(.5*n*(n-1)), (j-1)*n + i - .5*j*(j+1), 1)
  
}

# commutation matrix ====

#  DEF'N from Magnus & Neudecker (1979, p. 383, 3.1a)
comMat <- function(m, n=m){
  K <- matrix(0, nrow = n*m, ncol = n*m)
  for(i in 1:m){
    for(j in 1:n){
      H <- t(t(make_e_basis(i,m))) %*% t(make_e_basis(j,n))
      K <- K + kronecker(H, t(H))
    }
  }
  return(K)
}


# N Matrix ====

nMat <- function(n){
  return(.5 * (diag(n^2) + comMat(n)))
}  

# L Matrix ====

eliMat <- function(n, type = "s"){
  
  
  if (type == "s"){
    
    # case = "symmetric" (M vecA = vechA)
    #  DEF'N from Magnus & Neudecker (1980, p. 425, 3.1b)
    #  L in Magnus & Neudecker (1980)
    
    M <- matrix(0, nrow= .5*n*(n+1), ncol = n^2)
    
    for(i in 1:n){
      for(j in 1:n){
        if( i >= j){
          
          M <- M + kronecker(kronecker(make_u_basis(i,j,n),t(make_e_basis(j,n))),t(make_e_basis(i,n)))
          
        }
      }
    }
    
  } else if (type == "l"){
    
    # case = "strictly lower-triangular" (M vecA = vecpA)
    #  DEF'N from Neudecker (1983, p. 276, 3.1b)
    #  + L* in Neudecker (1983)
    
    M <- matrix(0, nrow= .5*n*(n-1), ncol = n^2)
    
    for(i in 1:n){
      for(j in 1:n){
        if( i > j){
          
          M <- M + kronecker(kronecker(make_v_basis(i,j,n),t(make_e_basis(j,n))),t(make_e_basis(i,n)))
          
        }
      }
    }
    
  } else if (type == "d"){
    
    # case = "diagonal" (M vecA = vecdA)
    #  DEF'N from Neudecker (1983, p. 279, 3.3a)
    #   + L** in Neudecker (1983)
    M <- matrix(0, nrow = n, ncol = n^2)
    
    for(i in 1:n){
      
      M <- M + kronecker(kronecker(make_e_basis(i,n), t(make_e_basis(i,n))), t(make_e_basis(i,n)))
      
    }
    
  } else {
    
    M <- NULL 
    
  }
  
  return(M)
  
}  

# duplication matrix ====

# DEF'N from Magnus & Neudecker (1980, p. 428, 3.2b)
dupMat <- function(n, type = "s"){
  
  if (type == "s"){
    
    # case = "symmetric" (D vechA = vecA)
    #   DEF'N from Magnus & Neudecker (1980, p. 428, 3.2b)
    # 
    
    Dt <- matrix(0, nrow = .5*n*(n+1), ncol= n^2)
    
    for(i in 1:n){
      for(j in 1:n){
        if( i >= j){
          Tm <- matrix(0,n,n); Tm[i,j] <- Tm[j,i] <- 1
          Dt <- Dt + make_u_basis(i,j,n) %*% t(vec(Tm))
        }
      }
    }
    
    return(t(Dt))
    
    # Neudecker (1983) Definition 3.2.a p. 278
  } else if (type == "l"){
    
    Dt <- matrix(0, nrow = .5*n*(n-1), ncol= n^2)
    
    for(i in 1:n){
      for(j in 1:n){
        if( i > j){
          T1 <- matrix(0,n,n); T1[i,j] <- 1
          T2 <- matrix(0,n,n); T2[j,i] <- 1
          Tm <- T1-T2
          Dt <- Dt + make_v_basis(i,j,n) %*% t(vec(Tm))
        }
      }
    }
    
    return(t(Dt))
    
    # Neudecker (1983) Definition 3.3.b Lemma 3.7 p. 280
  } else if( type == "d"){
    
    Dt <- eliMat(n, type = "d")
    return(t(Dt))
    
    
  }
  
}



makeCommutationMat <- function(n, m=n){
  
  #  DEF'N from Magnus & Neudecker (1979, p. 383, 3.1a)
  
  K <- matrix(0, nrow = n*m, ncol = n*m)
  
  for(i in 1:m){
    for(j in 1:n){
      
      H <- t(t(make_e_basis(i,m))) %*% t(make_e_basis(j,n))
      K <- K + kronecker(H, t(H))
      
    }
  }
  return(t(K))
}

genTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

matInv <- function(mat){
  
   minv <- genTryCatch(chol2inv(chol(mat)))
   
   if(is.null(minv$error)){
     
     return(minv$value)
     
   } else {
     
     return(solve(mat))
     
   }
   
}


