#' Summary information for a MIIV search object 
#'
#' @param object An object of class \code{miivs}
#' @param miivs.out A logical indicating whether the model-implied 
#'        instrumental variables found for \code{model} should be printed
#'        to the console. This is a temporary convenience function to 
#'        provide an editable, baseline format, for the \code{instruments}
#'        argument of \code{\link{miive}}.
#' @param eq.info A logical indicating whether equation-specific 
#'        information should be printed. Useful in models with a large number 
#'        of variables.
#' @param ... Optional arguments to summary, not used by user.
#' @export
summary.miivs <- function(object, miivs.out = FALSE, eq.info = FALSE,...){
  
  print(object)
  
  if (eq.info){
  
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
    
    invisible(lapply(object$eqns, function(eq){
      
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
          
        } else { 
          # do nothing
        }
      }
      cat("\n\n")
    }))
    
    cat("\n")
    cat("\n")
  }
  

  if(miivs.out){
    cat("instruments <- '\n")
      invisible(lapply(object$eqns,function(eq){
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