#' @method summary miive 
#' @export
summary.miive <- function(fit,eq.info = FALSE,...){
  
  
  fUp <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x
  }
  
  # Print 
  print(fit)
  
  cat("\n\nEquation Level Information:\n\n")
  
  w1 <- 20 # width of column 1
  w2 <- 20 # width of column 2

  if (eq.info){
    
    lbs <- c(
      "Equation Number",
      "Equation Type",
      "Equation Estimator",
      "Outcome Variable",
      "Explanatory Variable(s)",
      "Instrumental Variable(s)",
      "Categorical Variable(s)"
    )
    
    #et <- estimatesTable(fit)
    #eq <- fit$eqn[[1]]
    
    invisible(lapply(fit$eqn, function(eq){
      
      for(l in lbs){
        if(l == "Equation Number") {
          
          m1 <- paste0(l,": ")
          m2 <- paste0(eq$EQnum)
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
        } else if (l == "Equation Type") {
          
          m1 <- paste0(l,": ")
          m2 <- paste0(fUp(eq$EQmod))
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
        } else if (l == "Equation Estimator") {
          
          m1 <- paste0(l,": ")         
          m2 <- paste0(
            if (fit$estimator == "2SLS"){
              if (eq$categorical) {
                "MIIV-2SLS (PIV)"
              } else {
                "MIIV-2SLS"
              }
            } else if (fit$estimtor == "GMM"){
              "" # placeholder for GMM
            }
          )
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          cat("\n")
          
        } else if (l == "Outcome Variable") {
          
          msg <- paste0(l,": ",eq$DVlat)
          writeLines(strwrap(msg, indent = 0,exdent = 2, width = w1+w2))
          
        } else if (l == "Explanatory Variable(s)") {       
          
          msg <- paste0(l,": ",paste(eq$IVlat, collapse = ", "))
          writeLines(strwrap(msg, indent = 0,exdent = 2, width = w1+w2))

        } else if (l == "Instrumental Variable(s)") { 
          
          msg <- paste0(l,": ",paste(eq$MIIVs, collapse = ", "))
          writeLines(strwrap(msg, indent = 0,exdent = 2, width = w1+w2))
          
        } else if (l == "Categorical Variable(s)") { 
          
          cat.var <- paste0(
            if(eq$categorical){
              paste(intersect(
                c(eq$DVobs, eq$IVobs, eq$MIIVs), 
                fit$ordered
              ), collapse = ", ")            
            } else {
              "None"
            }
          )
          
          if (cat.var == "None") next
          
          msg <- paste0(l,": ",cat.var)
          writeLines(strwrap(msg, indent = 0,exdent = 2, width = w1+w2))
        } else { 
          # do nothing
        }
      }
      cat("\n\n")
    }))
    
  }
  

}