#' @export
print.miivs <- function(x,...){
  
  
  if (!x$eq.info){
    
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
  
     modeqns$'Composite Disturbance' <- NULL
    
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
    
  } else {
    
    fUp <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x
    }

    
    cat("\n\nMIIV SEARCH INFORMATION:\n\n")
    
    w1 <- 27 # width of column 1
    w2 <- 26 # width of column 2
    
    lbs <- c(
      "Equation Number",
      "Equation Type",
      "Outcome Variable",
      "Explanatory Variable(s)",
      "Composite Disturbance",
      "Instrumental Variable(s)"
    )
    
    invisible(lapply(x$eqns, function(eq){
      
      for(l in lbs){
        if(l == "Equation Number") {
          
          m1 <- paste0("   ",l,": ")
          m2 <- paste0(eq$EQnum)
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
        } else if (l == "Equation Type") {
          
          m1 <- paste0("   ",l,": ")
          m2 <- paste0(fUp(eq$EQmod))
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
        } else if (l == "Outcome Variable") {
          
          m1 <- paste0("   ",l,": ") 
          m2 <- if(eq$DVlat != eq$DVobs){
            paste0(eq$DVlat,"(",eq$DVobs,")")
          } else {
            paste0(eq$DVobs)
          }
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
        } else if (l == "Explanatory Variable(s)") {  
          
          m1 <- paste0("   ",l,": ") 
          m2 <- c(" ")
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
          m3 <- mapply(function (lat, obs){
            if(lat != obs){ 
              paste0(lat,"(",obs,")")
            } else { 
              paste0(obs)
            }
          }, lat = eq$IVlat, obs = eq$IVobs, SIMPLIFY = TRUE)
          
          msg <- paste0(m3, collapse = ", ")
          writeLines(strwrap(msg, indent = 5,exdent = 5, width = w1+w2))
          
        } else if (l == "Composite Disturbance") {
          m1 <- paste0("   ",l,": ") 
          m2 <- c(" ")
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
          m3  <- paste(eq$CD, collapse = ", ")   
          msg <- paste0(m3, collapse = ", ")
          writeLines(strwrap(msg, indent = 5,exdent = 5, width = w1+w2))
            

        } else if (l == "Instrumental Variable(s)") { 
          
          m1 <- paste0("   ",l,": ") 
          m2 <- c(" ")
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
          msg <- paste(eq$MIIVs, collapse = ", ")
          writeLines(strwrap(msg, indent = 5,exdent = 5, width = w1+w2))
          
        } else if (l == "Composite Disturbance") { 
          
          cat.var <- paste0(
       
              paste(intersect(
                c(eq$DVobs, eq$IVobs, eq$MIIVs), 
                object$ordered
              ), collapse = ", ")            
    
          )
        } else { 
          # do nothing
        }
      }
      cat("\n\n")
    }))
    
    cat("\n")
    cat("\n")
  }
  

  if(x$miivs.out){
    cat("instruments <- '\n")
      invisible(lapply(x$eqns,function(eq){
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