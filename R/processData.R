#' Process data.
#' 
#' @param data A data frame containing only those variables specified in the model syntax.
#' @param sample.cov Numeric matrix. A sample variance-covariance matrix. The rownames and colnames must contain the observed variable names.
#' @param sample.mean A sample mean vector.
#' @param sample.nobs Number of observations in the full data frame.
#' @param factor.vars A user-supplied vector of variable names for any categorical variables in \code{data}.
#' 
#' @details 
#' The internal function \code{processData} returns a list \code{g} 
#' which contains:
#' \describe{
#'  \item{\code{sample.cov} }{ 
#'    When the raw data contains only continous variables \code{sample.cov}
#'    is the sample covariance matrix for all variables indicated in the model
#'    syntax. \code{sample.cov} is calculated  as 
#'    \code{cov(data)*(nrow(data)-1)/nrow(data)}. When 
#'    the raw data contains categorical variables, and these 
#'    variables are also specified by the user in the 
#'    \code{factor.vars} argument of the \code{miive} function, 
#'    \code{sample.cov} is the sample covariance matrix for any
#'    continous variables included in the model syntax. If there are no 
#'    continuous variables in the \code{model} statement \code{sample.cov} is
#'    set to \code{NULL}.  If the raw data matrix contains missing data for any
#'    continous variables, the sample covariance matrix for the 
#'    continuous variables \code{sample.cov} is estimated with
#'    the \code{FIML} estimator. This covariance matrix is often
#'    referred to as the EM covariance matrix. When the data is a mix of
#'    continous and categorical variables, and missing
#'    data is present in the continous variables, the EM covariance
#'    matrix is estimated for the continous variables only.  
#'    In coefficient estimation this EM  sample covariance matrix is only 
#'    used if all variables in a given equation 
#'    (e.g. depdent, explanatory and instruments) are continous and no
#'    cross-equation restrictions are present. Otherwise the polychoric 
#'    correlation matrix is used. 
#'  }
#'  \item{\code{sample.mean} }{ 
#'    A vector of mean values for any continous variables in the dataset (in
#'    column order). Mean values for any categorical variable are given as
#'    \code{NA}. When missing values are present, the mean vector for any 
#'    continous variables is estimated using the options described above.
#'  }
#'  \item{\code{sample.nobs} }{ 
#'    The full number of rows contained in the user-supplied \code{data}, including
#'    those with missing data.
#'  }
#'  \item{\code{sample.polychoric} }{ 
#'    When the model syntax contains categorical variables the 
#'    \code{sample.polychoric} matrix is estimated using lavaan's
#'    \code{\link[lavaan]{lavCor}} function. If missing data is 
#'    present in the raw data matrix the polychoric correlations are 
#'    calculated using pairwise deletion. See the \code{\link[lavaan]{lavCor}} 
#'    documention for additional information. The \code{piv.opts} input 
#'    vector contains arguments that are passed to the 
#'    \code{estimator} and \code{se} arguments of \code{\link[lavaan]{lavCor}}.
#'    The default options for \code{piv.opts} are \code{estimator = "two.step"} 
#'    and \code{se = "standard"}.
#'  }
#'  \item{\code{asymptotic.cov} }{ 
#'    The asymptotic covariance matrix is estimated when either categorical
#'    variables are included in the model or when there is missing data. 
#'    Otherwise, \code{asymptotic.cov} is set to \code{NULL}. 
#'  }
#'  \item{\code{var.nobs} }{ 
#'    A numeric vector containing the number of non-missing observations 
#'    for each variable in \code{data}. If there are no missing observatuibs
#'    \code{var.nobs} is set to \code{NULL}.
#'  }
#'  \item{\code{var.categorical} }{ 
#'    A logical vector indicating which columns of \code{data} are categorical.
#'    If there are no categorical variables \code{var.categorical} is 
#'    set to \code{NULL}.
#'  }
#'  \item{\code{var.missing} }{ 
#'    A logical vector indicating which columns of \code{data} contain missing
#'    values.  If there are no variables with missing data
#'    \code{var.missing} is \code{NULL}.
#'  }
#' }
#' 
#' @examples
#' 
#'   
#'@keywords internal
processData <- function(data = data, 
                           sample.cov = sample.cov,
                           sample.mean = sample.mean, 
                           sample.nobs = sample.nobs, 
                           factor.vars = factor.vars){
  
  if (is.null(data)){
    
    if (!is.null(factor.vars)){
      stop(paste("miive: if categorical.vars are declared raw data is required."))
    }
    
    sample.sscp <- buildSSCP(sample.cov, sample.mean, sample.nobs)
    sample.polychoric <- NULL
    asymptotic.cov <- NULL
    var.nobs <- NULL
    var.categorical <- NULL
    var.missing <- NULL
    
  } else {
    
  
    piv.opts <- c(
      estimator = "two.step", 
      se = "standard"
    )

    missing.opts <- c(
      missing = "FIML", 
      estimator = "ML",         # se = std. implies m.v. normality
      se = "standard",          # or "robust.huber.white".
      information = "observed"  # or "expected"
    )
  
    sample.nobs     <- nrow(data)
    var.nobs        <- nrow(data) - colSums(is.na(data))
    var.missing     <- sapply(var.nobs, function(x) ifelse(x==sample.nobs, FALSE, TRUE))
    var.categorical <- vapply(data, is.factor, c(is.factor=FALSE))

    if( any(var.categorical) | !is.null(factor.vars) ){ 
      
      # if there are factors in data not given by the user
      
      if (is.null(factor.vars)){
        stop(paste("miive: undeclared factors in dataset.", 
                   "use the factor.vars argument to specify "))
      }
      
      #if (!identical(sort(colnames(data[,is.factor(data)])),sort(factor.vars))){
      #  stop(paste("miive: factor variables in data and factor.vars argument do not match."))
      #}
      
      fit <- lavaan::lavCor(
        data, 
        output = "fit", 
        missing = "pairwise",
        estimator = piv.opts["estimator"],
        se = piv.opts["se"],
        ordered = factor.vars 
      )
      
      # Polychoric correlation matrix. 
      sample.polychoric <- unclass(lavaan::inspect(fit, "cov.ov"))
      
      # Asymptotic covariance matrix of polychoric correlations. 
      asymptotic.cov  <- unclass(lavaan::inspect(fit, "vcov"))
      
      # Remove thresholds.
      asymptotic.cov  <- asymptotic.cov[
        1:(1/2*nrow(sample.polychoric)*(nrow(sample.polychoric))),
        1:(1/2*nrow(sample.polychoric)*(nrow(sample.polychoric)))
      ]

    } else {
      
      sample.polychoric <- NULL
      
    }
    
    continous.vars <- colnames(data)[!var.categorical]
    
    if (length(continous.vars > 1)){
      
      # Are there any missing observations
      if (any(var.nobs[continous.vars] < sample.nobs)){ # begin missing data
        
        var.cov <- outer(continous.vars, continous.vars, function(x, y) paste(x, "~~", y))
        model.saturated <- c(var.cov[lower.tri(var.cov, diag = TRUE)],paste(continous.vars, "~ 1"))
        
        fit <- lavaan::lavaan(
          model = model.saturated, 
          data = data[,continous.vars], 
          meanstructure = TRUE, 
          conditional.x = FALSE, 
          fixed.x = FALSE,
          missing = missing.opts["missing"], 
          estimator = missing.opts["estimator"], 
          se = missing.opts["se"], 
          information = missing.opts["information"]
        )
  
        # 2-step sample covariance matrix and mean vector.
        sample.cov  <- unclass(lavaan::lavInspect(fit, "cov.ov"))
        sample.mean[continous.vars] <- unclass(lavaan::lavInspect(fit, "mean.ov"))
       
        # Asymptotic covariance matrix of polychoric correlations. 
        asymptotic.cov  <- unclass(lavaan::inspect(fit, "vcov"))
        
        asymptotic.cov  <- asymptotic.cov[
          1:(1/2*nrow(sample.cov)*(nrow(sample.cov)+1)),
          1:(1/2*nrow(sample.cov)*(nrow(sample.cov)+1))
          ]
        
        sample.sscp <- buildSSCP(sample.cov, sample.mean, sample.nobs)
        
        # Remove covariances among means?
         # asymptotic.cov  <- asymptotic.cov[
         #   c(((1/2*nrow(sample.cov)*(nrow(sample.cov)+1))+1):nrow(asymptotic.cov),
         #     1:(1/2*nrow(sample.cov)*(nrow(sample.cov)+1))),
         #   c(((1/2*nrow(sample.cov)*(nrow(sample.cov)+1))+1):nrow(asymptotic.cov),
         #     1:(1/2*nrow(sample.cov)*(nrow(sample.cov)+1)))
         # ]
         # 
         # asymptotic.cov  <- rbind("1~1" = sample.nobs, cbind("1~1" = sample.nobs, asymptotic.cov))
   
      
      } else { # end missing data
        sample.cov  <- cov(data[,continous.vars])*(nrow(data[,continous.vars])-1)/nrow(data[,continous.vars])
        sample.mean <- colMeans(data[,continous.vars])
        sample.sscp <- buildSSCP(sample.cov, sample.mean, sample.nobs)
        if( is.null(factor.vars) ){ asymptotic.cov <- NULL }
      }
    } # end continuous only
  
  }
 
  # Prepare return list.
  g <- list(
    sample.cov  = sample.cov,
    sample.mean = sample.mean,
    sample.nobs = sample.nobs,
    sample.polychoric = sample.polychoric,
    sample.sscp = sample.sscp,
    asymptotic.cov = asymptotic.cov,
    var.nobs = var.nobs,
    var.categorical = var.categorical,
    var.missing = var.missing
  )

}


