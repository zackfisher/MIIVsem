#' Summary information for a MIIV estimation object
#'
#' @param object An object of class \code{miive}
#' @param eq.info A logical indicating whether equation-specific 
#'        information should be printed. Useful in models with 
#'        large numbers of variables.
#' @param restrict.tests A logical indicating whether two 
#'        test statistics for a large-sample wald test of 
#'        linaer hypotheses imposed on the MIIV-2SLS coefficient 
#'        matrix should be provided. The first statistic is 
#'        an approximate F and the second is Chi-square. 
#'        Assumptions and additional details for each test 
#'        are given by Greene (2000, p. 346-347) and Henningsen 
#'        and Hamman (2007).
#' @param rsquare A logical indicating whether R-square values for
#'        endogeneous variables are included in the output. Only
#'        available when \code{var.cov} is \code{TRUE}.
#' @param ... Optional arguments to summary, not used by user.
#' 
#' @references 
#' 
#' Greene, W. H. (2000). Econometric analysis. Upper Saddle River, N.J: 
#' Prentice Hall.
#' 
#' Henningsen, A., and Hamann, J.D. (2007). systemfit: A Package for 
#' Estimating Systems of Simultaneous Equations in R. Journal of Statistical 
#' Software 23(4), 1-40. 
#' 
#' @export
summary.miive <- function(object, eq.info = FALSE,
                          restrict.tests = FALSE,
                          rsquare = FALSE,...){
  
  fUp <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x
  }
  
  # Print 
  print(object)
  
  w1 <- 27 # width of column 1
  w2 <- 26 # width of column 2
  
  if (eq.info){
    
    cat("\n\nEQUATION INFORMATION:\n\n")
    
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
  
  if(rsquare & !is.null(object$v$rsquare) & object$var.cov){
    
    x   <- data.frame(
      "lhs" = names(object$v$rsquare),
      "op"  = "",
      "rhs"  = "",
      "est" = as.numeric(object$v$rsq),
      stringsAsFactors = FALSE
    )
     
 
    nd          <- 3L
    num.format  <- paste("%", max(8, nd + 5), ".", nd, "f", sep = "")
    char.format <- paste("%", max(8, nd + 5), "s", sep="") 
    
    y <- as.data.frame(
      lapply(x, function(x) {
        if(is.numeric(x)) {
          sprintf(num.format, x)   
        } else {
          x
        }
      }),
      stringsAsFactors = FALSE
    )
  
    y$op <- y$rhs<- NULL

    m <- as.matrix(
      format.data.frame(
        y, na.encode = TRUE, justify = "right"
      )
    )
  
    rownames(m) <- rep("", nrow(m))
    
      # rename some column names
    colnames(m)[ colnames(m) ==     "lhs" ] <- ""
    colnames(m)[ colnames(m) ==      "op" ] <- ""
    colnames(m)[ colnames(m) ==     "rhs" ] <- ""
    colnames(m)[ colnames(m) ==     "est" ] <- "Estimate"
    
    row.idx <- which(x$op == "")
    m[row.idx,1] <- makeName(x$lhs[row.idx])
    M <- m[row.idx,,drop=FALSE]
    colnames(M) <- colnames(m)
    rownames(M) <- rep("", NROW(M))
    cat("\n", "R-SQUARE", ":\n", sep = "")
    print(M, quote = FALSE)

  }
   
  if(restrict.tests & !is.null(object$r$R)){
    
    r.tests <- restrict.tests(object)
    
    cat("\n\nMIIV-2SLS LINEAR HYPOTHESIS TESTS:\n\n")
    
    r.rests <- rownames(object$r$R)
    for (i in 1:length(r.rests)){
      cat(paste0("   ",r.rests[i], "\n"))
    }

    cat("\n")
     
    m1 <- paste0("   ","Wald Test (Chi^2)",": ")
    m2 <- paste0(round(r.tests$chi.test,4))
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    # Degrees of Freedom
    m1 <- paste0("   ","Degrees of freedom",": ")
    m2 <- paste0(r.tests$chi.df)
    cat(sprintf("%-*s %*s\n", w1, m1, w2, m2));
    # P-value
    m1 <- paste0("   ","Pr(>Chi^2) ",": ")
    m2 <- paste0(round(r.tests$chi.p,4))
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