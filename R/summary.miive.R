#' @export
summary.miive <- function(object, eq.info = FALSE,...){
  
  fUp <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x
  }
  
  # Print 
  print(object)
  
  if (eq.info){
    
    cat("\n\nEQUATION INFORMATION:\n\n")
    
    w1 <- 27 # width of column 1
    w2 <- 26 # width of column 2
    
    lbs <- c(
      "Equation Number",
      "Equation Type",
      "Equation Estimator",
      "Outcome Variable",
      "Explanatory Variable(s)",
      "Instrumental Variable(s)",
      "Categorical Variable(s)"
    )
    
    invisible(lapply(object$eqn, function(eq){
      
      for(l in lbs){
        if(l == "Equation Number") {
          
          m1 <- paste0("   ",l,": ")
          m2 <- paste0(eq$EQnum)
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
        } else if (l == "Equation Type") {
          
          m1 <- paste0("   ",l,": ")
          m2 <- paste0(fUp(eq$EQmod))
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
        } else if (l == "Equation Estimator") {
          
          m1 <- paste0("   ",l,": ")       
          m2 <- paste0(
            if (object$estimator == "2SLS"){
              if (eq$categorical) {
                "MIIV-2SLS (PIV)"
              } else {
                "MIIV-2SLS"
              }
            } else if (object$estimtor == "GMM"){
              "" # placeholder for GMM
            }
          )
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

        } else if (l == "Instrumental Variable(s)") { 
          
          m1 <- paste0("   ",l,": ") 
          m2 <- c(" ")
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
          msg <- paste(eq$MIIVs, collapse = ", ")
          writeLines(strwrap(msg, indent = 5,exdent = 5, width = w1+w2))
          
        } else if (l == "Categorical Variable(s)") { 
          
          cat.var <- paste0(
            if(eq$categorical){
              paste(intersect(
                c(eq$DVobs, eq$IVobs, eq$MIIVs), 
                object$ordered
              ), collapse = ", ")            
            } else {
              "None"
            }
          )
          
          if (cat.var == "None") next
          
          m1 <- paste0("   ",l,": ") 
          m2 <- c(" ")
          cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
          
          writeLines(strwrap(cat.var, indent = 5,exdent = 5, width = w1+w2))
          
        } else { 
          # do nothing
        }
      }
      cat("\n\n")
    }))
    
  }
  

}