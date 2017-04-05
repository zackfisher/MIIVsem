#' Build the sum of squares and cross products matrix
#' 
#' Build sum of squares and crossproducts matrix (SSCP).
#' From the means, covariances, and n's you can recover the
#' raw sum-of-squares and products matrix for all the variables. 
#' Say the matrix of all the variables is X, with mean vector bar(x), 
#' and covariance matrix S, based on sample-size n. Then the SSCP 
#' matrix is X'X = (n - 1)S + n bar(x) bar(x)'. You then need to 
#' add the row/column for the constant, which is just n in the 1, 1 
#' position and n bar(x) elsewhere.
#' 
#' @param sample.cov Numeric matrix. A sample variance-covariance matrix. The rownames and colnames must contain the observed variable names.
#' @param sample.mean A sample mean vector.
#' @param sample.nobs Number of observations in the full data frame.
#' 
#' @details 
#' \describe{
#'  \item{\code{sample.cov} }{ 
#'    When the raw data contains only continous variables \code{sample.cov}
#'    is the sample covariance matrix for all variables indicated in the model
#'    syntax calculated  as \code{cov(data)*(nrow(data)-1)/nrow(data)}. When 
#'    the raw data contains categorical variables, and these 
#'    variables are also specified by the user in the 
#'    \code{factor.vars} argument of the \code{miive} function, 
#'    \code{sample.cov} is the sample covariance matrix for any
#'    continous variables included in the model syntax. If there are no 
#'    continuous variables in the \code{model} statement \code{sample.cov} is
#'    set to \code{NULL}.  If the raw data matrix contains missing data for any
#'    continous variables the sample covariance matrix for the 
#'    continuous variables in the data are estimated using the saturated model 
#'    and the \code{missing = FIML} option. When the data is a mix of
#'    continous and categorical variables and missing
#'    data is present in the continous variables, the saturated model is estimated
#'    using the continous variables only.  In coefficient estimation this EM 
#'    sample covariance matrix is only used if all variables in a given 
#'    equation (e.g. depdent, explanatory and instruments) are continous and no
#'    cross-equation restrictions are present. Otherwise the polychoric correlation 
#'    matrix is used. NEED MORE INFORMATION HERE.
#'  }
#'  \item{\code{sample.mean} }{ 
#'    A vector of mean values for any continous variables in the dataset (in
#'    column order). Mean values for any categorical variable are given as
#'    \code{NA}. When missing values are present the mean vector for any continous
#'    variables is estimated using the options described above.
#'  }
#'  \item{\code{sample.nobs} }{ 
#'    The number of rows contained in the user-supplied \code{data}.
#'  }
#' }
#'@keywords internal
buildSSCP <- function(sample.cov, sample.mean, sample.nobs){
  
  if (is.null(sample.cov)){
    
    res <- NULL
    
  } else {
    
    res <- matrix(NA, length(sample.mean)+1, length(sample.mean)+1,
                  dimnames = list(c("1", names(sample.mean)),
                                  c("1", names(sample.mean))))
  
    res[1,1]   <- sample.nobs
    res[-1,1]  <- res[1,-1] <- sample.mean*sample.nobs
    res[-1,-1] <- (sample.cov + sample.mean %*% t(sample.mean))*sample.nobs
  }
  
  return(res)
  
}