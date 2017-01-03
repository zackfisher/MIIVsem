#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using model-implied instrumental variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by as.data.frame to data frame.
#' @param instruments A user-supplied list of valid MIIVs for each equation. See Example 2 below. 
#' @param estimator Options \code{"2SLS"} or \code{"GMM"} for estimating the model parameters. Default is \code{"2SLS"}.
#' @param control .
#'
#' @return model
#' @return dat
#' @return modeeqns
#' 
#' @details 
#' \itemize{
#' \item{\code{instruments}} {Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.  To utilize this option you must first define a list of instruments using the syntax displayed below. After the list is defined, set the \code{instruments} argument equal to the name of the list of MIIVs. Note, \code{instruments} are specified for an equation, and not for a specific endogenous variable.}
#' }
#' 
#' @references 
#' 
#' Bollen, K. A. 1996.	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 
#' Bollen,	K. A. 2001.	Two-stage	Least	Squares	and	Latent	Variable	Models:	
#' Simultaneous	Estimation	and	Robustness	to	Misspecifications.
#' In	R.	Cudeck,	S.	Du	Toit,	and	D.	Sorbom	(Eds.),	Structural	
#' Equation	Modeling:	Present	and	Future,	A	Festschrift	in	Honor	of	Karl	
#' Joreskog	(pp. 119-138).	Lincoln,	IL: Scientific	Software.
#' 	
#'
#' @examples
#' 
#' # Example 1
#' 
#'  bollen1989a_model <- '
#'
#'    Eta1 =~ y1 + y2  + y3  + y4  
#'    Eta2 =~ y5 + y6  + y7  + y8    
#'    Xi1  =~ x1 + x2 + x3 
#'
#'    Eta1 ~ Xi1  
#'    Eta2 ~ Xi1 
#'    Eta2 ~ Eta1 
#'
#'    y1   ~~ y5
#'    y2   ~~ y4
#'    y2   ~~ y6
#'    y3   ~~ y7
#'    y4   ~~ y8
#'    y6   ~~ y8 
#'  '
#'  
#'   miive(model = bollen1989a_model, data = bollen1989a)
#'  
#'  
#' # Example 2
#' 
#'   my_instruments <- ' 
#'    y1 ~ x2 + x3                            
#'    y5 ~ y2 + y3 + y4 + x2                
#'    y2 ~ y3 + y7 + y8 + x2           
#'    y3 ~ y2 + y4 + y6 + y8        
#'    y4 ~ y3 + y6           
#'    y6 ~ y3 + y4 + y7 + x2            
#'    y7 ~ y2 + y4 + y6 + y8       
#'    y8 ~ y2 + y3 + y7 + x2          
#'    x2 ~ y1 + y5 + y2 + y3 + y4 + y6
#'    x3 ~ y1 + y5 + y2 + y3 + y4 + y6
#'  '
#'  
#' miive(model = bollen1989a_model, data = bollen1989a, 
#'       instruments = my_instruments)
#'  
#'  
#' # Example 3
#'  bollen1989a_model_r <- '
#'
#'    Eta1 =~ y1 + l2*y2  + l3*y3  + l4*y4  
#'    Eta2 =~ y5 + l2*y6  + l3*y7  + l4*y8    
#'    Xi1  =~ x1 + x2 + 0.5*x3 
#'
#'    Eta1 ~ Xi1  
#'    Eta2 ~ Xi1 
#'    Eta2 ~ Eta1 
#'
#'    y1   ~~ y5
#'    y2   ~~ y4
#'    y2   ~~ y6
#'    y3   ~~ y7
#'    y4   ~~ y8
#'    y6   ~~ y8
#'    
#'  '
#'  
#'  miive(model = bollen1989a_model_r, data = bollen1989a)
#'  
#'  
#' @export
miive <- function(model = model, data = NULL, instruments = NULL, 
                  estimator = "2SLS", control = NULL){
  
  #-------------------------------------------------------#  
  # Check class of model.
  #-------------------------------------------------------#
  if ( "miivs" == class(model) ){ mod <- model; d <- mod$eqns; } 
  if ( "miivs" != class(model) ){ mod <- miivs(model); d <- mod$eqns; } 
  
  #-------------------------------------------------------# 
  # generateFormulas
  #-------------------------------------------------------# 
  #   input:  (1) miivs equation list 'd' (miivs(foo)$eqns)
  #                returned from miivs search function.
  #           (2) instruments argument, null by default,
  #               if no instruments have been supplied.
  #          
  #  return:  (1) updated 'd' (miivs(foo)$eqns)
  #
  #  details: 'd' updated with three new fields:
  #           (1) EqFormula is two-sided formula for equation i.
  #           (2) MIIVsFormula is a one-sided formula
  #           characterizing the instruments for equation i.
  #           (3) MIIVsUsed is a vector of MIIVs to be used
  #           in estimation for equation i.
  #
  #           User specified instruments can be provided
  #           using the 'instruments' argument. The 'instruments'
  #           object should be a character of similar form to the
  #           model syntax except that only the "~" operator is 
  #           needed.  For the dependent variables, Z1 and Z2 we
  #           could specify the instruments M1, M2,and M3 as follows:
  #
  #           instruments <- '
  #             Z1 ~ M1 + M2 + M3
  #             Z2 ~ M1 + M2 + M3  
  #           '
  #           
  #           The function contains checks to determine if (1)
  #           the user-supplied dependent variables exist in the
  #           set of estimating equations and (2) if the instruments
  #           provided by the user are valid MIIVs.
  #-------------------------------------------------------#
  d     <- generateFormulas(d, instruments)
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Prepare data.  (REVISIT)
  #   
  #  To-do: Prepare for missing data.
  #         Prepare for covariance matrices
  #-------------------------------------------------------#
  mf <- prepareRawData(data)
  #-------------------------------------------------------#

  #-------------------------------------------------------#  
  # Prepare lists for terms and matrices.
  #-------------------------------------------------------#
  dvLabels    <- unlist(lapply(d,"[[","DVobs"))
  numEq       <- length(d)      # number of equations
  numCoef_i   <- numeric(numEq) # vector of number of coefficients in each eq.
  coefNames_i <- list() # any names of coefficients of each equation
  obsNames_i  <- list() # any names of observations of each equation
  terms_zy_i  <- list() # list of terms for Z and Y objects in each equation
  terms_v_i   <- list() # list of terms for V objects in each equation
  Y_i  <- list()        # list of dependent variable vectors in each eq.
  Z_i  <- list()        # list of explanatory variable matrices in each eq.
  V_i  <- list()        # list of instrumental variable matrices in each eq.
  
  for(i in 1:numEq) { 
    
    mf$formula <- NULL; evalMF <- NULL;
    mf$formula   <- d[[i]]$EqFormula
    evalMF       <- eval(mf)
    
    terms_zy_i[[i]] <- attr( evalMF, "terms" )
    Y_i[[i]]        <- model.extract(evalMF, "response" )
    Z_i[[i]]        <- model.matrix(terms_zy_i[[i]], evalMF)
    
    obsNames_i[[i]]  <- rownames(Z_i[[i]])
    
    numCoef_i[i]     <- ncol(Z_i[[i]])
    
    coefNames_i[[i]] <- paste0(d[[i]]$DVobs,"_", colnames(Z_i[[i]]))
    
    mf$formula <- NULL; evalMF <- NULL;
    mf$formula   <- d[[i]]$MIIVsFormula
    evalMF       <- eval(mf)
    
    terms_v_i[[i]] <- attr(evalMF, "terms")
    V_i[[i]] <- model.matrix(terms_v_i[[i]], evalMF)
  }

  numCoef      <- sum(numCoef_i)
  coefNames    <- unlist(coefNames_i)

  #-------------------------------------------------------#  
  # Process any missing values.
  #-------------------------------------------------------#
  noNA        <- processMissing(Y_i, Z_i, V_i)
  numObsPerEq <- nrow(noNA)
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # For now remove any missing values.
  #-------------------------------------------------------#
  for(i in 1:numEq) {
    Y_i[[i]]   <- Y_i[[i]][noNA[,i]]
    attrAssign <- attributes(Z_i[[i]])$assign
    Z_i[[i]]   <- Z_i[[i]][noNA[, i], , drop = FALSE]
    attributes(Z_i[[i]])$assign <- attrAssign
    V_i[[i]] <- V_i[[i]][noNA[, i], , drop = FALSE ]
  }
  
  validObsNames_i <- lapply(Y_i,names)
  numObsEq_i      <- unlist(lapply(Y_i,length))
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Calculate stage 1 fitted regression matrices using 
  # raw data with the missing values removed.  
  #-------------------------------------------------------#
  Zhat_i <- calcStage1Fitted(Z_i, V_i)
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Build Restriction Matrices.
  # returns NULL if there are no restrictions,
  # otherwise returns a list containing the 'R' matrix
  # and 'q' vector, as well as a vector 'cons' of 
  # the constrained coefficients.
  #-------------------------------------------------------#  
  restrictions <- buildRestrictMat(d)
  if (is.null(restrictions)){ R <- NULL; q <- NULL
  } else { R <- restrictions$R; q <- restrictions$q }
  #-------------------------------------------------------#   
  
  #-------------------------------------------------------#  
  # Calculate degrees of freedom.
  #-------------------------------------------------------#  
  df_list       <- calcDegreesOfFreedom(d, restrictions, numObsEq_i, numCoef_i)
  df_i          <- df_list$df
  numCoefRes    <- df_list$numCoefRes   # nCoefLiAll
  numCoefRes_i  <- df_list$numCoefRes_i # nCoefLiEq
  #-------------------------------------------------------# 
  
  #-------------------------------------------------------#  
  # Calculate stage 2 regression.
  #-------------------------------------------------------#
  b <- calcStage2Coefs(Y_i, Zhat_i, V_i, R, q, coefNames)
  #-------------------------------------------------------#
  
  #-------------------------------------------------------#  
  # Calculate the residuals.
  #-------------------------------------------------------# 
  resids_i <- calcResiduals(b, Y_i, Z_i, numCoef_i)
  #-------------------------------------------------------#   
    
  #-------------------------------------------------------#  
  # Calculate the residual covariance matrix.
  #-------------------------------------------------------#  
  Omega <- calcResidCovMat(resids_i,numEq, noNA = noNA) 
  #-------------------------------------------------------# 

  #-------------------------------------------------------#  
  # Calculate coefficient covariance matrix.
  #-------------------------------------------------------#  
  coefCov <- calcCoefCovMat(Zhat_i, Omega, R = R, q = q, numEq, noNA, coefNames)
  #-------------------------------------------------------#  
  
  #-------------------------------------------------------#  
  # Fitted.values
  #-------------------------------------------------------#  
  fitted.values <- bdiag(Z_i) %*% b  
  #-------------------------------------------------------#  
  
  #-------------------------------------------------------#  
  # Organize results
  #-------------------------------------------------------#
  b_i <- split(b, rep(1:length(numCoef_i), numCoef_i))
  results <- list()
  
  for(i in 1:numEq) {
    results$eq[[i]]          <- list()
    results$eq[[i]]$eqnNum   <- i         # equation number
    results$eq[[i]]$eqnLabel <- NA        # eqnLabels[[i]]
    results$eq[[i]]$method   <- estimator
    
    # Residuals
    results$eq[[i]]$residuals        <- resids_i[[i]]
    names(results$eq[[i]]$residuals) <- validObsNames_i[[i]]
    
    # Coefficients
    results$eq[[i]]$coefficients         <- b_i[[i]]
    #names(results$eq[[i]]$coefficients)  <- coefNames_i[[i]]
    
    # Coefficient Covariance Matrices
    results$eq[[i]]$coefCov  <- as.matrix(
      coefCov[(1+sum(numCoef_i[1:i])-numCoef_i[i]):(sum(numCoef_i[1:i])),
              (1+sum(numCoef_i[1:i])-numCoef_i[i]):(sum(numCoef_i[1:i]))] 
    )
    colnames(results$eq[[i]]$coefCov)    <- coefNames_i[[i]]
    rownames(results$eq[[i]]$coefCov)    <- coefNames_i[[i]]
    
    # Fitted Values
    results$eq[[i]]$fitted.values <- rep(NA, numObsPerEq)
    results$eq[[i]]$fitted.values[noNA[,i]] <-
      fitted.values[(1+sum(numObsEq_i[1:i])-numObsEq_i[i]):(sum(numObsEq_i[1:i]))]
    
    names(results$eq[[i]]$fitted.values) <- obsNames_i[[i]]
    
    # Terms
    results$eq[[i]]$terms <- terms_zy_i[[i]]
    
    # Rank
    results$eq[[i]]$rank  <- numCoefRes_i[i]
    
    # Total # of coef. in system
    results$eq[[i]]$nCoef.sys <- numCoef
    
    # Total # of rest. coef. in system
    results$eq[[i]]$rank.sys  <-  numCoefRes    # nCoefLiAll
    
    # rank = number of linear independent coefficients of the entire system
    results$eq[[i]]$df.residual  <- df_i[i]    
    results$eq[[ i ]]$df.residual.sys  <- sum(numObsEq_i) - numCoefRes
    results$eq[[i]]$inst      <- d[[i]]$MIIVsUsed
    results$eq[[i]]$termsInst <- terms_v_i[[i]]
    class(results$eq[[i]]) <- "miive.equation"
  }

  # System coefficients
  results$coefficients <- as.numeric(drop(b))
  names(results$coefficients) <- coefNames
  
  # Full Coefficient Covariance Matrix
  results$coefCov <- as.matrix(coefCov)
  dimnames(results$coefCov) <- list(coefNames, coefNames)
  
  # Residual Covariance Matrix 
  results$residCovEst <- as.matrix(Omega)
  dimnames(results$residCovEst) <- list(dvLabels, dvLabels)

  # Residual Covarance Matrix
  results$residCov <- function(resids_i, numEq, noNA, oneSigma = FALSE,
                               centered = FALSE, diagOnly = FALSE)
  dimnames(results$residCov) <- list(dvLabels, dvLabels)
  
  results$method  <- estimator
  results$rank    <- numCoefRes
  # rank = total number of linear independent coefficients of all equations
  results$df.residual <- sum(numObsEq_i) - numCoefRes
  # degrees of freedom of the whole system
  results$restrict.matrix <- R
  results$restrict.rhs    <- q
  results$control <- control
  class(results)  <- "miive"
  
  results
}
