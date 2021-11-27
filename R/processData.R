#' Process data.
#'@keywords internal
processData <- function(data = data, 
                        sample.cov = sample.cov,
                        sample.mean = sample.mean, 
                        sample.nobs = sample.nobs, 
                        ordered = ordered,
                        missing = missing,
                        se = se,
                        pt = pt){
  
  # if the user supplied a covariance matrix.
  if (is.null(data)){
    
    if (is.null(sample.mean) ){
      sample.mean <- rep(0L, ncol(sample.cov))
      names(sample.mean) <- colnames(sample.cov)
    }
    
    if (!is.null(ordered)){
      stop(paste("miive: if categorical.vars are declared raw data is required."))
    }
    
    sample.sscp        <- buildSSCP(sample.cov, sample.mean, sample.nobs)
    sample.polychoric  <- NULL
    asymptotic.cov     <- NULL
    asymptotic.cov.sat <- NULL
    var.nobs           <- NULL
    var.categorical    <- NULL
    var.missing        <- NULL
    var.exogenous      <- NULL
  
  # the user supplied raw data
  } else {
    
    asymptotic.cov     <- NULL
    asymptotic.cov.sat <- NULL
    
    # Data-level characteristics
    sample.nobs     <- nrow(data)
    
    # Variable level characteristics
    var.categorical <- vapply(data, is.factor, c(is.factor=FALSE))
    var.nobs        <- nrow(data) - colSums(is.na(data))
    var.missing     <- sapply(var.nobs, function(x) {
      ifelse(x==sample.nobs, FALSE, TRUE)
    })
 
    # Exogenous variables in the dataset
    var.exogenous <- colnames(data) %in%
      pt[ pt$exo == 1L & pt$op  == "~~" & pt$lhs == pt$rhs, "rhs" ]
    
    names(var.exogenous) <- colnames(data)
    
    #-------------------------------------------------------# 
    # Ordered factors present in data
    #-------------------------------------------------------# 
    
    if( any(var.categorical) | !is.null(ordered) ){ 
      
      # there are factors in the user-upplied data
      if (any(var.categorical)) { 
        
        ## if there are undeclared factors throw an error
        if (length(setdiff(colnames(data)[var.categorical],ordered)) > 0){
          und.factors <- setdiff(colnames(data)[var.categorical],ordered)
          stop(paste0(
            "miive: the following undeclared factors,",
            "were found in data: ", 
            paste(und.factors,collapse = ", "), 
            ". Use the ordered argument to specify",
            "categorical variables."
          ))
        }
      }
      
      # For now we don't do anything about ov.exogenous variables.
      ov.names.x <- NULL
      
      if (any(var.categorical) & missing != "listwise") { 
        stop(paste0(
          "miive: missing = ", missing, 
          " not supported in the presence of",
          " categorical variables."
        ))
      }
      
      
      fit <- lavaan::lavCor(
        data, 
        output = "fit", 
        missing = "listwise",
        estimator = "two.step",
        se = "standard",
        ov.names.x = ov.names.x,
        ordered = setdiff(ordered, ov.names.x)
      )
      
      sample.sscp <- NULL
      sample.nobs <- nrow(data)
 
      # Polychoric correlation matrix. 
      sample.polychoric <- unclass(lavaan::inspect(fit, "cor.ov"))
      
      # is polycor faster
      # corS<-matrix(NA,12,12)
      # for (i in 1:12){    
      #   for (j in i:12) {
      #     if(i==j){
      #       corS[i,j] <- 1
      #     } else {
      #       corS[i,j] <- polycor::polychor(data[,i],data[,j],ML=F)
      #       corS[j,i] <- corS[i,j]
      #     }
      #   }
      # }
      
      # Asymptotic covariance matrix of polychoric correlations. 
      asymptotic.cov  <- unclass(lavaan::inspect(fit, "vcov"))
        
      ordered.varnames <- apply(
        t(utils::combn(colnames(sample.polychoric), 2)), 1, function(x){
          paste0(x[1], "~~", x[2])
      })
        
        # Reorder asymptotic covariance matrix.
      asymptotic.cov  <- asymptotic.cov[
        ordered.varnames, 
        ordered.varnames
      ]
     
    } else { # there are no categorical variables
      
      sample.polychoric <- NULL
      
    }
    
    continuous.vars <- colnames(data)[!var.categorical]
    
    if (length(continuous.vars) > 1){
      
      # Are there any missing observations
      # any.miss <- any(var.nobs[continuous.vars] < sample.nobs)
      
      if (missing == "twostage"){ # begin missing data
        
        var.cov <- outer(
          continuous.vars, continuous.vars, 
          function(x, y) {
            paste(x, "~~", y)
          }
        )
        
        saturated.model <- c(
          var.cov[lower.tri(var.cov, diag = TRUE)],
          paste(continuous.vars, "~ 1")
        )
        
        saturated.fit <- lavaan::lavaan(
          model = saturated.model, 
          data = data[,continuous.vars], 
          meanstructure = TRUE, 
          conditional.x = FALSE, 
          fixed.x = FALSE,
          missing = "FIML", 
          estimator = "ML", 
          se = "robust.huber.white", 
          information = "observed"
        )
        
        #lavInspect(fit, "sampstat.h1")
        sample.cov  <- unclass(lavaan::inspect(saturated.fit, "cov.ov"))
        sample.mean <- unclass(lavaan::lavInspect(saturated.fit, "mean.ov"))
        sample.nobs <- lavaan::lavInspect(saturated.fit, "nobs") 
        sample.sscp <- buildSSCP(sample.cov, sample.mean, sample.nobs)
        asymptotic.cov.sat <- unclass(lavaan::vcov(saturated.fit))

      } else { # end missing data
        
        sample.cov  <- stats::cov(data[,continuous.vars])*
          (nrow(data[,continuous.vars])-1) / 
          nrow(data[,continuous.vars])
        
        sample.mean <- colMeans(data[,continuous.vars])
        sample.nobs <- nrow(data)
        sample.sscp <- buildSSCP(sample.cov, sample.mean, sample.nobs)
      
      }
    } 
  
  }
 
  # Prepare return list.
  g <- list(
    sample.cov  = sample.cov,
    sample.mean = sample.mean,
    sample.nobs = sample.nobs,
    sample.polychoric = sample.polychoric,
    sample.sscp = sample.sscp,
    asymptotic.cov = asymptotic.cov,
    asymptotic.cov.sat = asymptotic.cov.sat,
    var.nobs = var.nobs,
    var.categorical = var.categorical,
    var.missing = var.missing,
    var.exogenous = var.exogenous
  )

}


