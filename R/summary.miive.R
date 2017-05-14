#' @export
summary.miive <- function(object, eq.info = FALSE,
                          restrict.tests = FALSE,...){
  
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
  

  if(restrict.tests){
    
    r.tests <- restrict.tests(object)
    
    cat("\n\nLINEAR HYPOTHESIS TESTS:\n\n")
    
    r.rests <- rownames(object$r$R)
    for (i in 1:length(r.rests)){
      cat(paste0("   ",r.rests[i], "\n"))
    }

    cat("\n")
     
    m1 <- paste0("   ","Wald Test (Chi^2)",": ")
    m2 <- paste0(round(r.tests$wald.test,4))
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    # Degrees of Freedom
    m1 <- paste0("   ","Degrees of freedom",": ")
    m2 <- paste0(r.tests$wald.df)
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    # P-value
    m1 <- paste0("   ","Pr(>Chi^2) ",": ")
    m2 <- paste0(round(r.tests$wald.p,4))
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    
    
    cat("\n")
     
    m1 <- paste0("   ","Wald Test (F)",": ")
    m2 <- paste0(round(r.tests$f.test,4))
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    # Degrees of Freedom
    m1 <- paste0("   ","Degrees of freedom",": ")
    m2 <- paste0(paste(c(r.tests$f.df1,r.tests$f.df1), collapse=", "))
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    # P-value
    m1 <- paste0("   ","Pr(>F) ",": ")
    m2 <- paste0(round(r.tests$f.p,4))
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    
  }
}