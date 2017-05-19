#' Print method for a MIIV estimation object
#' 
#' @param x a miive object
#' @param ... Optional arguments to print, not used by user.
#' 
#' @export
print.miive <- function(x,...){
  
  ## an ugly simplification of the truly
  ## impressive lavan::lav_print function

  # What esimtators were used?
  eq.estimators <- unlist(lapply(x$eqn, function(eq){
    if (x$estimator == "2SLS"){
      if (eq$categorical){
        "MIIV-2SLS (PIV)"
      } else {
        "MIIV-2SLS"
      }
    }
  }))
  
  disp.estimator <- if (all(eq.estimators == "MIIV-2SLS (PIV)")){
    "MIIV-2SLS (PIV)"
  } else if (all(eq.estimators == "MIIV-2SLS")){
    "MIIV-2SLS"
  } else {
    "Mixed"
  }


  # MIIVsem version number
  cat(paste0("MIIVsem (", utils::packageVersion("MIIVsem"),") results"), "\n\n")
  
  w1 <- 40 # width of column 1
  w2 <- 36 # width of column 2
  
  head.txt  <- do.call("rbind",
    list(
      c("Number of observations", x$sample.nobs),
      c("Number of equations", length(x$eqn)),
      c("Estimator", disp.estimator),
      c("Standard Errors", ifelse(
        any(x$se %in% c("boot", "bootstrap")), "bootstrap", x$se
      )),
      c("Missing", x$missing)
    )
  )
  for(i in 1:nrow(head.txt)){
    cat(sprintf("%-*s %*s\n", w1, head.txt[i,1], w2, head.txt[i, 2]));
  }
    
  cat("\n")
  
  if (x$se  %in% c("boot", "bootstrap")){
    
    if (x$boot.ci == "norm"){ boot.ci <- "Normal"}     else if
       (x$boot.ci == "perc"){ boot.ci <- "Percentile"} else if
       (x$boot.ci == "basic"){boot.ci <- "Basic"}      else if
       (x$boot.ci == "bca"){  boot.ci <- "Adjusted Percentile"} else 
       {                      boot.ci <- x$boot.ci            }
    
    boot.txt  <- do.call("rbind",
      list(
        c("Bootstrap reps requested", x$bootstrap),
        c("Bootstrap reps successful", x$bootstrap.true),
        c("Bootstrap intervals", boot.ci)
      )
    )
    for(i in 1:nrow(boot.txt)){
      cat(sprintf("%-*s %*s\n", w1, boot.txt[i,1], w2, boot.txt[i, 2]));
    }
    cat("\n")
  } 
  

  
  if(length(x$eqn.unid) > 0){
    cat(paste0("***", "Insufficient MIIVs to estimate the ",
      paste(lapply(x$eqn.unid, "[[","DVlat"),collapse=", "),
    " equations.***"), "\n")
  }
    
    
  sections <-  c("MEASUREMENT MODEL", 
                 "STRUCTURAL MODEL", 
                 "INTERCEPTS",
                 "VARIANCES",
                 "COVARIANCES")
  
  x    <- estimatesTable(x, sarg = TRUE)
  
  ## remove duplicate Sargan test info from 
  # regression equations to avoid confusion
  dup.row.nums <- row.names(
    x[x$op == "~",][which(duplicated(x[x$op == "~",]$lhs)),]
  )
  
  if (length(dup.row.nums) > 0) {
    x[dup.row.nums, c("sarg", "sarg.df", "sarg.p")] <- NA
  }

  x$eq <- NULL
  
  nd          <- 3L
  num.format  <- paste("%", max(8, nd + 5), ".", nd, "f", sep = "")
  char.format <- paste("%", max(8, nd + 5), "s", sep="") 
  
  cat("\nParameter Estimates:\n\n")
  
  # round to 3 digits after the decimal point
  y <- as.data.frame(
    lapply(x, function(x) {
      if(is.numeric(x)) {
        sprintf(num.format, x)   
      } else {
        x
      }
    }),
    stringsAsFactors = FALSE)
  
  # fix degress of freedom for saragan test
  suppressWarnings(
    y$sarg.df    <- as.numeric(as.character(y$sarg.df))
  )
  sarg.format  <- paste("%", max(4), ".", 0, "f", sep = "")
  y$sarg.df    <- sprintf(sarg.format, y$sarg.df)

  
  y$op <- y$rhs<- NULL

  m <- as.matrix(
    format.data.frame(
      y, na.encode = TRUE, justify = "right"
    )
  )
  
  rownames(m) <- rep("", nrow(m))
  
  if(!is.null(x$sarg)) {

    sarg.idx <- which(is.na(x$sarg))
    
    if(length(sarg.idx) > 0L) {
      m[sarg.idx, "sarg"] <- ""
      
      if(!is.null(x$sarg.df)) {
        m[sarg.idx, "sarg.df"] <- ""
      }
      if(!is.null(x$sarg.p)) {
        m[sarg.idx, "sarg.p"] <- ""
      }
    }
  }
  
    if(!is.null(x$lower)) {

      boot.idx <- which(is.na(x$lower)) 
    
      if(length(boot.idx) > 0L) {
      
        if(!is.null(x$lower)) {
          m[boot.idx, "lower"] <- ""
        }
        
        if(!is.null(x$lower)) {
          m[boot.idx, "upper"] <- ""
        }
    }
  }
  
  if(!is.null(x$se)) {
    
    se.idx <- which(is.na(x$se))
    
    if(length(se.idx) > 0L) {
      m[se.idx, "se"] <- ""
      
      if(!is.null(x$z)) {
        m[se.idx, "z"] <- ""
      }
      if(!is.null(x$pvalue)) {
        m[se.idx, "pvalue"] <- ""
      }
      if(!is.null(x$upper)) {
        m[se.idx, "upper"] <- ""
      }
      if(!is.null(x$lower)) {
        m[se.idx, "lower"] <- ""
      }
    }
    
    se.idx <- which(x$se == 0)
    
    if(length(se.idx) > 0L) {
      m[se.idx, "se"] <- ""
      if(!is.null(x$z)) {
        m[se.idx, "z"] <- ""
      }
      if(!is.null(x$pvalue)) {
        m[se.idx, "pvalue"] <- ""
      }
      if(!is.null(x$upper)) {
        m[se.idx, "upper"] <- ""
      }
      if(!is.null(x$lower)) {
        m[se.idx, "lower"] <- ""
      }
    }
  }
  
  # rename some column names
  colnames(m)[ colnames(m) ==     "lhs" ] <- ""
  colnames(m)[ colnames(m) ==      "op" ] <- ""
  colnames(m)[ colnames(m) ==     "rhs" ] <- ""
  colnames(m)[ colnames(m) ==     "est" ] <- "Estimate"
  colnames(m)[ colnames(m) ==      "se" ] <- "Std.Err"
  colnames(m)[ colnames(m) ==       "z" ] <- "z-value"
  colnames(m)[ colnames(m) ==  "pvalue" ] <- "P(>|z|)"
  colnames(m)[ colnames(m) ==   "lower" ] <- "Lower"
  colnames(m)[ colnames(m) ==   "upper" ] <- "Upper"
  colnames(m)[ colnames(m) ==    "sarg" ] <- "Sargan"
  colnames(m)[ colnames(m) == "sarg.df" ] <- "df"
  colnames(m)[ colnames(m) ==  "sarg.p" ] <- "P(Chi)"
  
  #colnames(m) <- sprintf(char.format,  colnames(m))
  colnames(m)[grepl("df", colnames(m))] <- sprintf(
    paste("%", max(4), "s", sep=""),  
    colnames(m)[grepl("df", colnames(m))]
  )
  colnames(m)[!grepl("df", colnames(m))] <- sprintf(
    char.format,  
    colnames(m)[!grepl("df", colnames(m))]
  )

  ht <- TRUE
  
  for(s in sections) {
    if(s == "MEASUREMENT MODEL") {
      row.idx <- which(x$op == "=~" & x$rhs != "1")
      if(length(row.idx) == 0L) next
      m[row.idx,1] <- makeName(x$rhs[row.idx])
    } else if(s == "STRUCTURAL MODEL") {
      row.idx <- which( x$op == "~" & x$rhs != "1")
      if(length(row.idx) == 0L) next
      m[row.idx,1] <- makeName(x$rhs[row.idx])
    } else if(s == "COVARIANCES") {
      row.idx <- which(x$op == "~~" & x$lhs != x$rhs)
      if(length(row.idx) == 0L) next
      m[row.idx,1] <- makeName(x$rhs[row.idx])
    } else if(s == "INTERCEPTS") {
      row.idx <- which(x$op  == "~1")
      if(length(row.idx) == 0L) next
      m[row.idx,1] <- makeName(x$lhs[row.idx])
    } else if(s == "VARIANCES") {
      row.idx <- which(x$op == "~~" & x$lhs == x$rhs)
      if(length(row.idx) == 0L) next
      m[row.idx,1] <- makeName(x$rhs[row.idx])
    } else {
      row.idx <- integer(0L)
    }
    
    if(s %in% c("MEASUREMENT MODEL",
                "STRUCTURAL MODEL", 
                "COVARIANCES")) {
      
      nel <- length(row.idx)
      
      M <- matrix("", nrow = nel*2, ncol = ncol(m))
      
      colnames(M) <- colnames(m)
      rownames(M) <- rep("", NROW(M))
     
      LHS <- paste(x$lhs[row.idx], x$op[row.idx])
      lhs.idx <- seq(1, nel*2L, 2L)
      rhs.idx <- seq(1, nel*2L, 2L) + 1L
      PREFIX  <- rep("", length(LHS))
      M[lhs.idx, 1] <- sprintf("%1s%-15s", PREFIX, LHS)
      M[rhs.idx,  ] <- m[row.idx,]
      # avoid duplicated LHS labels
      if(nel > 1L) {
        del.idx <- integer(0)
        old.lhs <- ""
        for(i in 1:nel) {
          if(LHS[i] == old.lhs) {
            del.idx <- c(del.idx, lhs.idx[i])
          }
          old.lhs <- LHS[i]
        }
        if(length(del.idx) > 0L) {
          M <- M[-del.idx,,drop=FALSE]
        }
      }
      
      if(s %in% c("MEASUREMENT MODEL", "STRUCTURAL MODEL")) {
        if(ht){
          cat("\n", "STRUCTURAL COEFFICIENTS", ":\n", sep = "")
          ht <- FALSE
        } else {
          #cat("\n\n", sep = "")
          colnames(M)<- rep("", length(colnames(M)))
        }
      } else {
        cat("\n", s, ":\n", sep = "")
      }
      
    
      if(s %in% "COVARIANCES"){
        colnames(M)[ grepl("Sargan",  colnames(M))] <- ""
        colnames(M)[ grepl("df",      colnames(M))] <- ""
        colnames(M)[ grepl("Chi",    colnames(M)) ] <- ""
        print(M, quote = FALSE)
      } else {
        print(M, quote = FALSE)
      }
     
        
    } else {
      M <- m[row.idx,,drop=FALSE]
      colnames(M) <- colnames(m)
      rownames(M) <- rep("", NROW(M))
      cat("\n", s, ":\n", sep = "")
      if(s %in% c("INTERCEPTS", "VARIANCES")){
        
        M[,grepl("Sargan",  colnames(M))] <- ""
        M[,grepl("df",      colnames(M))] <- ""
        M[,grepl("Chi",     colnames(M))] <- ""
        
        colnames(M)[ grepl("Sargan",  colnames(M))] <- ""
        colnames(M)[ grepl("df",      colnames(M))] <- ""
        colnames(M)[ grepl("Chi",     colnames(M))] <- ""
        print(M, quote = FALSE)
      } else {
        print(M, quote = FALSE)
      }
    }
  } 

}