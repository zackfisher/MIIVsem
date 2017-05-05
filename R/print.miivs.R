#' @export
print.miivs <- function(x,...){
  
  z <- x$eqns
  
  for (i in 1:length(z)){
    LHS <- paste(z[[i]]$DVobs, collapse = ", ")
    RHS <- paste(z[[i]]$IVobs, collapse = ", ")
    Instruments <- paste(z[[i]]$MIIVs, collapse = ", ")
    Disturbance <- paste(z[[i]]$CDist, collapse = ", ", sep="")
    modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
    colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
    if (i == 1) {modeqns <- modtemp }
    if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
  }
  
  if (!x$composite.disturbance){
    modeqns$'Composite Disturbance' <- NULL
  }
  
  cat("Model Equation Information \n")
  cat("\n")
  print(
    modeqns, 
    quote = FALSE, 
    right = FALSE, 
    row.names = FALSE, 
    print.gap=1
  )
  cat("\n")
  cat("\n")
  
  if(x$miivs.out){
    cat("instruments <- '\n")
      invisible(lapply(z,function(eq){
        cat(
          paste0(
            "   ",
            eq$DVobs,
            " ~ ",
            paste(eq$MIIVs, collapse = " + "),
            "\n"
          )
        )
      }))
    cat("'")
  }

}