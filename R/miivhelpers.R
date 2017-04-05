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
  r <- unlist(sapply(eqns,"[[", c("coefficients")))
  z <- cbind(do.call(rbind, strsplit(names(r), "~")), r)
  z <- z[which(z[,2]!="1"),] # cuts down on nrow
  
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
  
  return(mod)
}
