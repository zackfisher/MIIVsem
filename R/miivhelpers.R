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

revVech <- function (x, nrow = NULL){
  dim(x) <- NULL
  if (is.null(nrow)) 
    nrow <- (-1 + sqrt(1 + 8 * length(x)))/2
  output <- matrix(0, nrow, nrow)
  output[lower.tri(output, diag = TRUE)] <- x
  hold <- output
  hold[upper.tri(hold, diag = TRUE)] <- 0
  output <- output + t(hold)
  return(output)
}