#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate structural equation models using model-implied instrumental 
#' variables (MIIVs).
#'
#' @param model A model specified using lavaan model syntax or a
#'        \code{\link{miivs}} object returned by the \code{\link{miivs}}
#'        function. See Details for more information about permissible
#'        operators and example model syntax.
#' @param data A data frame, list or environment or an object coercible 
#'        by \code{as.data.frame} to data frame. The most common 
#'        application is to supply a data.frame.
#' @param instruments This allows user to specify the instruments for 
#'        each equation. See Details and the \code{miivs.out} argument of 
#'        \code{\link{summary.miivs}} for more information on the correct 
#'        input format. External (auxiliary) instruments can be 
#'        supplied, however, the \code{miiv.check} argument must
#'        be set to \code{FALSE}. In the typical application, the program 
#'        will choose the MIIVs for each equation based on the model 
#'        specification. To view the model implied instruments after 
#'        estimation see the \code{eq.info} argument of 
#'         \code{\link{summary.miive}}.
#' @param sample.cov Numeric matrix. A sample variance-covariance 
#'        matrix. The rownames and colnames attributes must contain all
#'        the observed variable names indicated in the model syntax.
#' @param sample.mean A sample mean vector. If \code{sample.cov} is provided
#'        and the \code{sample.mean} argument is \code{NULL}, intercepts
#'        for all endogenous variables will not be estimated. 
#' @param sample.nobs Number of observations in the full data frame.
#' @param sample.cov.rescale If \code{TRUE}, the sample covariance matrix
#'        provided by the user is internally rescaled by multiplying 
#'        it with a factor (N-1)/N.
#' @param estimator Options \code{"2SLS"} or \code{"GMM"} for estimating the
#'        model parameters. Default is \code{"2SLS"}. Currently, only 
#'        \code{2SLS} is supported.
#' @param est.only If \code{TRUE}, only the coefficients are returned.
#' @param var.cov If \code{TRUE}, variance and covariance parameters are
#'        estimated.
#' @param se If "standard", asymptotic standard errors are
#'        computed. If \code{"bootstrap"} (or \code{"boot"}), bootstrap 
#'        standard errors are computed.
#' @param bootstrap Number of bootstrap draws, if bootstrapping is used. The
#'        default is \code{1000}
#' @param boot.ci Method for calculating bootstrap confidence intervals. 
#'        Options are normal approximation (\code{"norm"}), basic bootstrap 
#'        interval (\code{"basic"}), percentile interval (\code{"perc"}), 
#'        and adjusted bootstrap percentile (\code{"bca"}). The default is 
#'        normal approximation. See \code{\link[boot]{boot.ci}} for more 
#'        information. 
#' @param missing Default is \code{"listwise"} however, a maximum likelihood 
#'        related missing data method called \code{"twostage"} is 
#'        also available. See Details below  on \code{missing} for more 
#'        information.
#' @param miiv.check Default is \code{TRUE}. \code{miiv.check} provides a 
#'        check to determine whether user-upplied instruments are implied
#'        by the model specification (e.g. valid MIIVs). When auxiliary or 
#'        external instruments are provided \code{miiv.check} should be 
#'        set to \code{FALSE}. 
#' @param ordered A vector of variable names to be treated as ordered factors
#'        in generating the polychoric correlation matrix and subsequent PIV
#'        estimates. See details on \code{ordered} below for more information.
#' @param sarg.adjust Adjusment methods used to adjust the p-values associated
#'        with the Sargan test due to multiple comparisons. Defaults is 
#'        \code{none}. For options see \code{\link[stats]{p.adjust}}.
#' @param overid.degree A numeric value indicating the degree of 
#'        overidentification to be used in estimation. 
#' @param overid.method The method by which excess MIIVs should
#'        be pruned to satisfy the \code{overid.degree}. Options include
#'        random (\code{minimum.eigen}) or stepwise R2.
#'        (\code{stepwise.R2}).The default is \code{stepwise.R2}.
#' @param information The type of information to be used in the two-step 
#'        missing data approach. Options are "observed" or "expected".
#' @param twostage.se The type of standard error to be used in the two-step 
#'        missing data approach. Options are "standard" or "robust.huber.white".
#' @param auxiliary Names of any auxiliary variables to be used in the two-step 
#'        missing data approach. 
#' @details 
#' 
#' \itemize{
#' \item{\code{model}} {
#' 
#'   The following model syntax operators are currently supported: =~,
#'   ~, ~~ and *. See below for details on default behavior, descriptions 
#'   of how to specify the scaling indicator in latent variable models, 
#'   and how to impose equality constraints on the parameter estimates. 
#'   
#'   \strong{Example using Syntax Operators}
#'   
#'   In the model below, 'L1 =~ Z1 + Z2 + Z3'  indicates the latent variable L1
#'   is measured by 3 indicators, Z1, Z2, and Z3. Likewise, L2 is measured by 3
#'   indicators, Z4, Z5, and Z6. The statement 'L1 ~ L2' specifies latent
#'   variable L1 is regressed on latent variable L2. 'Z1 ~~ Z2' indicates the
#'   error of Z2 is allowed to covary with the error of Z3. The label LA3
#'   appended to Z3 and Z6 in the measurement model constrains the factor
#'   loadings for Z3 and Z6 to equality. For additional details on constraints
#'   see Equality Constraints  and Parameter Restrictions.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + Z2 + LA3*Z3
#'      L2 =~ Z4 + Z5 + LA3*Z6
#'      L1  ~ L2
#'      Z2 ~~ Z3
#'   '}
#'  
#'   \strong{Scaling Indicators}
#'   
#'   Following the \pkg{lavaan} model syntax, latent variables are defined using
#'   the =~ operator.  For first order factors, the scaling indicator chosen is
#'   the first observed variable on the RHS of an equation. For the model  below
#'   \code{Z1} would be chosen as the scaling indicator for \code{L1} and
#'   \code{Z4} would be chosen as the scaling indicator for \code{L2}.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + Z2 + Z3
#'      L2 =~ Z4 + Z5 + Z6
#'   '}
#'   
#'   \strong{Equality Constraints and Parameter Restrictions}
#'   
#'   Within- and across-equation equality constraints on the factor loading and
#'   regression coefficients can be imposed directly in the model syntax. To
#'   specify equality constraints between different parameters equivalent labels
#'   should be prepended to the variable name using the * operator. For example,
#'   we could constrain the factor loadings for the two non-scaling indicators
#'   of \code{L1} to equality using the following  model syntax.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + LA2*Z2 + LA2*Z3
#'      L2 =~ Z4 + Z5 + Z6
#'   '}
#'   
#'   Researchers also can constrain the factor loading and regression 
#'   coefficients to specific numeric values in a similar fashion. Below we
#'   constrain the regression coefficient  of \code{L1} on \code{L2} to
#'   \code{1}.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + Z2 + Z3
#'      L2 =~ Z4 + Z5 + Z6
#'      L3 =~ Z7 + Z8 + Z9
#'      L1  ~ 1*L2 + L3
#'   '}
#'
#'   \strong{Higher-order Factor Models}
#'   
#'   For example, in the model below, the  scaling indicator for the 
#'   higher-order factor \code{H1} is taken to be \code{Z1}, the scaling 
#'   indicator that would have been assigned to the first lower-order factor
#'   \code{L1}. The intercepts for lower-order latent variables are set to zero,
#'   by default
#'   
#'   \preformatted{model <- '
#'      H1 =~ L1 + L2 + L3
#'      L1 =~ Z1 + Z2 + Z3
#'      L2 =~ Z4 + Z5 + Z6
#'      L3 =~ Z7 + Z8 + Z9
#'   '}
#'   
#'   \strong{Model Defaults}
#'   
#'   In addition to those relationships specified in the model syntax 
#'   \pkg{MIIVsem} will automatically include the intercepts of any observed or
#'   latent endogenous variable. The intercepts for any scaling indicators and
#'   lower-order latent variables are set to zero by default. Covariances among
#'   exogenous latent and observed  variables are included when \code{var.cov =
#'   TRUE}. Where appropriate the covariances of the errors of latent and
#'   observed dependent variables are automatically included in the model
#'   specification. These defaults correspond to those used by \pkg{lavaan} and
#'   \code{auto = TRUE}, except that endogenous latent variable intercepts are
#'   estimated by default, and the intercepts of scaling indicators are fixed to
#'   zero.
#'   
#'   
#'   \strong{Invalid Specifications}
#'   
#'   Certain model specifications are not currently supported.  For example, the
#'   scaling indicator of a latent variable is not permitted to cross-load on
#'   another latent variable. In the model below \code{Z1}, the scaling
#'   indicator for L1, cross-loads on the latent variable \code{L2}. Executing a
#'   search on the model below will result in the warning: \emph{miivs: scaling
#'   indicators with a factor complexity greater than 1 are not currently
#'   supported}.
#'   
#'   \preformatted{model <- '
#'     L1 =~ Z1 + Z2 + Z3
#'     L2 =~ Z4 + Z5 + Z6 + Z1
#'   '}
#'   
#'   In addition, \pkg{MIIVsem} does not currently support relations
#'   where the scaling indicator of a latent variable is also the 
#'   dependent variable in a regression equation.  The
#'   model below would not be valid under the current algorithm.
#'   
#'   \preformatted{model <- '
#'     L1 =~ Z1 + Z2 + Z3
#'     Z1  ~ Z4
#'     Z4  ~ Z5 + Z6
#'   '}
#'   }
#'   
#'   \item{\code{instruments}} {
#' 
#'   To utilize this option you must first define a list of instruments using
#'   the syntax displayed below. Here, the dependent variable for each equation
#'   is listed on the LHS of the ~ operator. In the case of latent variable
#'   equations, the dependent variable is the scaling indicator associated with
#'   that variable. The instruments are then given on the RHS, separated by +
#'   signs. The instrument syntax is then encloses in single quotes. For example,
#'   
#'   \preformatted{customIVs <- '
#'      y1 ~ z1 + z2 + z3
#'      y2 ~ z4 + z5
#'   '
#'   }
#'     
#'   After this list is defined, set the \code{instruments} argument equal to 
#'   the name of the list of instruments (e.g. \code{customIVs}). Note, that
#'   \code{instruments} are specified for an equation, and not for a specific
#'   endogenous variable. If only a subset of dependent variables are listed in
#'   the instruments argument, only those  equations listed will be estimated. 
#'   If external or auxiliary instruments (instruments not otherwise included in
#'   the model) are included the \code{miiv.check} argument should be set to
#'   \code{FALSE}.
#'   }
#'   
#'   \item{\code{sample.cov}} {
#'   
#'   The user may provide a sample covariance matrix in lieu of raw data. The
#'   rownames and colnames must contain the observed variable names indicated in
#'   the model syntax. If \code{sample.cov} is not \code{NULL} the user must
#'   also supply a vector of sample means (\code{sample.mean}), and the number
#'   of sample observations (\code{sample.nobs}) from which the means and
#'   covariances were calculated.  If no vector of sample means is provided 
#'   intercepts will not be estimated. \pkg{MIIVsem} does not support bootstrap 
#'   standard errors or polychoric instrumental variable estimtation when the
#'   sample moments, rather than raw data, are used as input.
#'   }
#'   
#'   \item{\code{sample.mean}} {
#'   
#'   A vector of length corresponding to the row and column dimensions
#'   of the \code{sample.cov} matrix.  The names of \code{sample.mean} 
#'   must match those in the \code{sample.cov}.  If the user supplies a 
#'   covariance matrix but no vector of sample means intercepts will not 
#'   be estimated.
#'   }
#'   
#'   \item{\code{sample.cov.rescale}} {
#'   
#'   Default is \code{TRUE}.  If the sample covariance matrix provided 
#'   by the user should be internally rescaled by multiplying it with a 
#'   factor (N-1)/N.
#'   }
#'   
#'   \item{\code{estimator}} {
#'   
#'   The default estimator is \code{2SLS}. For equations with continuous 
#'   variables only and no restrictions the estimates are identical to those
#'   described in Bollen (1996, 2001). If restrictions are present a restricted
#'   MIIV-2SLS estimator is implemented using methods similar to those described
#'   by Greene (2003) but adapted for moment based estimation. 2SLS coefficients
#'   and overidentifcation tests are constructed using the sample moments for 
#'   increased computational efficiency.
#'   
#'   If an equation contains ordered categorical variables, declared in the
#'   \code{ordered} argument, the PIV estimator described by Bollen and
#'   Maydeu-Olivares (2007) is implemented. The PIV estimator does not currently
#'   support exogenous observed predictors of endogenous categorical variables. 
#'   See details of the \code{ordered} argument for more information about the
#'   PIV estimator.
#'   }
#'   
#'  \item{\code{se}} {
#'   When \code{se} is set to \code{"boot"} or \code{"bootstrap"} standard 
#'   errors are computed using a nonparametric bootstrap assuming an independent
#'   random sample. If \code{var.cov = TRUE} nonceonvergence may occur and any
#'   datasets with impproper solutions will be recorded as such and discarded.
#'   Bootstrapping is implemented using the \pkg{boot} by resampling the
#'   observations in \code{data} and refitting the model with the resampled
#'   data. The number of bootstrap replications is set using the
#'   \code{bootstrap} argument, the default is \code{1000}. Here, the standard
#'   errors are based on the standard deviation of successful bootstrap
#'   replications.   Note, the Sargan test statistic is calculated from the
#'   original sample and is not a bootstrap-based estimate. When \code{se} is
#'   set to \code{"standard"} standard errors for the MIIV-2SLS coefficients are
#'   calculated using analytic expressions. For equations with categorical
#'   endogenous variables, the asymptotic distribution of the coefficients is
#'   obtained via a first order expansion where the matrix of partial
#'   derivatives is evaluated at the sample polychoric correlations. For some
#'   details on these standard errors see Bollen & Maydeu-Olivares (2007, p.
#'   315). If \code{var.cov = TRUE} only point estimates for the variance and 
#'   covariance estimates are calculated.  To obtain standard errors for the
#'   variance and covariance parameters we recommend setting \code{se =
#'   "bootstrap"}. Analytic standard errors for the variance covariance
#'   parameters accounting for the first stage estimation have been derived and
#'   will be available in future releases.
#'   }
#'   
#'  \item{\code{missing}} {
#'   There are two ways to handle missing data in \pkg{MIIVsem}. First, missing 
#'   data may be handled by listwise deletion (\code{missing = "listwise"}), In
#'   this case any row of data containing missing observations is excluded from
#'   the analysis and the sample moments are adjusted accordingly. Estimation
#'   then proceeds normally. The second option for handling missing data is
#'   through a two-stage procedures \code{missing = "twostage"} where consistent
#'   estimates of the saturated population means and covariances are obtained in
#'   the first stage. These quantities are often referred to as the "EM means"
#'   and "EM covariance matrix." In the second stage the saturated estimates are
#'   used to calculate the MIIV-2SLS structural coefficients. Bootstrap standard
#'   errors are recommended but will be  computationally burdensome due to the
#'   cost of calculating the EM-based moments at each bootstrap replication.
#'   }
#'  \item{\code{ordered}} {  
#'   For equations containing ordered categorical variables MIIV-2SLS
#'   coefficients are estimated using the approach outlined in Bollen
#'   & Maydeu-Olivares (2007). The asymptotic distribution of the 
#'   these coefficients is obtained via a first order expansion where 
#'   the matrix of partial derivatives is evaluated at the sample 
#'   polychoric correlations. For some details on these
#'   standard errors see Bollen & Maydeu-Olivares (2007, p. 315). If 
#'   \code{var.cov = TRUE} only point estimates for the variance and
#'   covariance estimates are calculated using the \code{ULS} estimator
#'   in \pkg{lavaan}. To obtain standard errors for the variance and 
#'   covariance parameters we recommend the bootstrap approach. 
#'   Analytic standard errors for the variance covariance parameters 
#'   in the presence of endogenous categorical variables 
#'   will be available in future releases. Currently \pkg{MIIVsem}
#'   does not support exogenous variables in equations with categorical
#'   endogenous variables.
#'   }
#' }
#' 
#' \strong{Sargan's Test of Overidentification}  
#' 
#' An essential ingredient in the MIIV-2SLS approach is the application of 
#' overidentification tests when a given model specification leads to an excess
#' of instruments. Empirically, overidentification tests are used to evalulate
#' the assumption of orthogonality between the instruments and equation 
#' residuals. Rejection of the null hypothesis implies a deficit in the logic
#' leading to the instrument selection. In the context of MIIV-2SLS this is the
#' model specification itself.  By default, \pkg{MIIVsem} provides Sargan's
#' overidentification test (Sargan, 1958) for each overidentified equation in
#' the system. When cross-equation restrictions or missing data are present the 
#' properties of the test are not known. When the system contains many equations
#' the \code{sarg.adjust} option provides methods to adjust the p-values
#' associated with the Sargan test due to multiple comparisons. Defaults is 
#' \code{none}. For other options see \code{\link[stats]{p.adjust}}.
#' 
#' @references 
#' 
#' Bollen, K. A. (1996).	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 
#' Bollen, K. A. (2001).	Two-stage	Least	Squares	and	Latent	Variable	
#' Models: Simultaneous	Estimation	and	Robustness	to	Misspecifications.
#' In	R.	Cudeck,	S.	Du	Toit,	and	D.	Sorbom	(Eds.),	Structural	
#' Equation	Modeling:	Present	and	Future,	A	Festschrift	in	Honor	of	Karl	
#' Joreskog	(pp. 119-138).	Lincoln,	IL: Scientific	Software.
#' 	
#' Bollen, K. A., & Maydeu-Olivares, A. (2007). A Polychoric Instrumental 
#' Variable (PIV) Estimator for Structural Equation Models with Categorical 
#' Variables. \emph{Psychometrika}, 72(3), 309.
#' 
#' Freedman, D. (1984). On Bootstrapping Two-Stage Least-Squares Estimates 
#' in Stationary Linear Models. \emph{The Annals of Statistics}, 
#' 12(3), 827–842. 
#' 
#' Greene, W. H. (2000). Econometric analysis. Upper Saddle River, N.J: 
#' Prentice Hall.
#' 
#' Hayashi, F. (2000). Econometrics. Princeton, NJ: Princeton University 
#' Press
#' 
#' Sargan, J. D. (1958). The Estimation of Economic Relationships 
#' using Instrumental Variables. Econometrica, 26(3), 393–415. 
#' 
#' Savalei, V. (2010). Expected versus Observed Information in SEM with 
#' Incomplete Normal and Nonnormal Data. \emph{Psychological Methods}, 
#' 15(4), 352–367. 
#' 
#' Savalei, V., & Falk, C. F. (2014). Robust Two-Stage Approach 
#' Outperforms Robust Full Information Maximum Likelihood With 
#' Incomplete Nonnormal Data. \emph{Structural Equation Modeling: 
#' A Multidisciplinary Journal}, 21(2), 280–302. 
#'
#' 
#' @seealso \link{MIIVsem}{miivs}
#' 
#' @keywords MIIV-2SLS MIIV PIV 2sls tsls instrument SEM two-stage least-squares
#'  
#' @export
miive <- function(model = model, 
                  data = NULL,  
                  instruments = NULL,
                  sample.cov = NULL, 
                  sample.mean = NULL, 
                  sample.nobs = NULL, 
                  sample.cov.rescale = TRUE, 
                  estimator = "2SLS", 
                  se = "standard", 
                  bootstrap = 1000L, 
                  boot.ci = "norm",
                  missing = "listwise",
                  est.only = FALSE, 
                  var.cov = FALSE, 
                  miiv.check = TRUE, 
                  ordered = NULL,
                  sarg.adjust = "none",
                  overid.degree = NULL,
                  overid.method = "stepwise.R2",
                  information = "observed",
                  twostage.se = "standard",
                  auxiliary = NULL
                  ){
  
  #-------------------------------------------------------#
  # In the current release disable "twostage" missing
  #-------------------------------------------------------#
  
  
  #-------------------------------------------------------# 
  # A few basic sanity checks for user-supplied covariance 
  # matrices and mean vectors. Also, rescale cov.matrix if
  # cov.rescale set to TRUE.
  #-------------------------------------------------------# 
  if(is.null(data)){
    
    if (!is.null(ordered)){
      stop(paste(
        "miive: raw data required when declaring factor variables.")
      )
    }
    
    if(!is.null(sample.mean)){
      
      if (!is.null(sample.mean) & !is.vector(sample.mean)){
        stop(paste(
          "miive: sample.mean must be a vector.")
        )
      }
      
      if (!all.equal( names(sample.mean), colnames(sample.cov), 
                    check.attributes = FALSE )){
        stop(paste(
          "miive: names of sample.mean vector",
          "and sample.cov matrix must match.")
        )
      }
      
    }
    
    if(sample.cov.rescale & !is.null(sample.cov)){
      sample.cov <- sample.cov * (sample.nobs-1)/sample.nobs
    }
    
  }
  
  #-------------------------------------------------------#  
  # Check class of model.
  #-------------------------------------------------------#
  if ( "miivs" == class(model) ){ 
    
    d  <- model$eqns 
    pt <- model$pt
    
  } else { 
    
    res <- miivs(model)
    d   <- res$eqns
    pt  <- res$pt
    
  } 
  
  #-------------------------------------------------------# 
  # parseInstrumentSyntax
  #  - primarily for user-supplied instruments
  #-------------------------------------------------------#
  d  <- parseInstrumentSyntax(d, instruments, miiv.check)
  
  #-------------------------------------------------------# 
  # remove equations where there are not sufficient MIIVs
  #-------------------------------------------------------#
  underid   <- unlist(lapply(d, function(eq) {
    length(eq$MIIVs) < length(eq$IVobs)
  }))
  d.un <- d[underid]; d   <- d[!underid]
  
  #-------------------------------------------------------# 
  # Remove variables from data that are not in model 
  # syntax and preserve the original column ordering.
  #-------------------------------------------------------# 
  if(!is.null(data)){
    
    obs.vars <- unique(unlist(lapply(d,"[", c("DVobs", "IVobs", "MIIVs"))))
    
    if(!is.null(auxiliary)){
      obs.vars <- unique(c(obs.vars, auxiliary))
    }
    
    if (any(!obs.vars %in% colnames(data))){
     stop(paste(
       "miive: model syntax contains variables not in data.")
     )
    }
    
    data <- data[,colnames(data) %in% obs.vars]
    data <- as.data.frame(data)
    
    # convert any variables listed in 
    # ordered to categorical
    data[,ordered]  <- lapply(
      data[,ordered, drop = FALSE], ordered
    )
  }
  
  #-------------------------------------------------------# 
  # If missing == "listwise"
  #-------------------------------------------------------# 
  if (!is.null(data) & missing == "listwise"){
    
    data <- data[stats::complete.cases(data),]
    
  } 

  #-------------------------------------------------------# 
  # Process data. See documentation of processData. 
  #-------------------------------------------------------# 
  g <- processData(data, 
                   sample.cov, 
                   sample.mean, 
                   sample.nobs, 
                   ordered, 
                   missing,
                   se,
                   pt,
                   information,
                   twostage.se,
                   auxiliary)
  
  #-------------------------------------------------------# 
  # Instrument pruning takes place here. 
  #-------------------------------------------------------# 
  if (!is.null(overid.degree)){
    d <- pruneExcessMIIVs(d, 
                          overid.degree,  
                          overid.method, 
                          data = data, 
                          sample.cov = g$sample.cov, 
                          sample.polychoric = g$sample.polychoric,
                          sample.mean = g$sample.mean, 
                          sample.nobs = g$sample.nobs,
                          cat.vars    = g$var.categorical)
  }
  
  #-------------------------------------------------------# 
  # Add some fields to d and check for any problematic
  # cases prior to estimation.  
  #-------------------------------------------------------# 
  d <- lapply(d, function (eq) {
    
    eq$missing <- ifelse(
      any(g$var.nobs[c(eq$DVobs, eq$IVobs, eq$MIIVs)] < g$sample.nobs), 
      TRUE, 
      FALSE
    )
    
    eq$categorical <- ifelse(
      any(g$var.categorical[c(eq$DVobs, eq$IVobs, eq$MIIVs)]), 
      TRUE, 
      FALSE
    )
    
    eq$exogenous   <- ifelse(
      any(g$var.exogenous[c(eq$DVobs, eq$IVobs, eq$MIIVs)]), 
      TRUE, 
      FALSE
    )
    
    # For now throw an error if a categorical equation
    # contains an exogenous variable. 
    if (eq$categorical & eq$exogenous){
      
      stop(paste("miive: exogenous variables in",
                 "equations containing categorical",
                 "variables (including MIIVs)",
                 "are not currently supported."))
    }
    # if (eq$missing & eq$restricted){
    #   stop(paste("miive: Restrictions on coefficients in",
    #              "equations with missing values",
    #              "(including on the MIIVs)",
    #              "is not currently supported."))
    # }
    eq
  })
  
  #-------------------------------------------------------#  
  # Build Restriction Matrices.
  #-------------------------------------------------------#  
  r <- buildRestrictMat(d)
  
  d <- lapply(d, function (eq) {

    eq$restricted  <- ifelse(
      r$eq.restricted[eq$DVobs],
      TRUE,
      FALSE
    )
    # throw an error if a categorical equation 
    # includes restrictions
    if ((eq$categorical & eq$restricted) & se == "standard"){
      stop(paste("miive: When SE = 'standard', restrictions",
                 "on coefficients in equations containing",
                 "categorical variables (including MIIVs)",
                 "are not allowed. Remove restrictions or",
                 "use se = 'bootstrap'."))
    }
    eq
  })
  
  #-------------------------------------------------------#
  # Generate results object.
  #-------------------------------------------------------#
  results <- switch(
    estimator,
      "2SLS" = miive.2sls(d, g, r, est.only, se, missing, var.cov, sarg.adjust),
      # "GMM"  = miive.gmm(d, d.un, g, r, est.only, se), # Not implemented
      # In other cases, raise an error
      stop(
        paste(
          "Invalid estimator:", 
          estimator, 
          "Valid estimators are: 2SLS"
        )
      )
  )
  
  #-------------------------------------------------------#
  # Point estimates for variance and covariance params
  #-------------------------------------------------------#
  if (var.cov){
    
    if(length(d.un) > 0){
      stop(paste("MIIVsem: variance covariance esimtation not",
                 "allowed in the presence of underidentified",
                 "equations."))
    }
    
  
    v <- estVarCovarCoef(data = data, 
                         g = g, 
                         eqns = results$eqn, 
                         pt = pt, 
                         ordered = ordered,
                         missing = missing)
    
  } else {
    
    v <- NULL
    
  }
  
  
  #-------------------------------------------------------#
  # Boostrap and substitute closed form SEs with boot SEs
  #-------------------------------------------------------# 
  results$boot <- NULL
  
  if(se == "boot" | se == "bootstrap"){
    
    if(is.null(data)){
      stop(paste("miive: raw data required for bootstrap SEs."))
    }
    
    # Do this outside of the bootstrap loop

    boot.results <- boot::boot(data, function(origData, indices){
      
      bsample <- origData[indices,]
      
      g.b <- processData(
        data = bsample, 
        sample.cov = NULL, 
        sample.mean = NULL, 
        sample.nobs = NULL, 
        ordered = ordered, 
        missing = missing,
        se = se,
        pt = pt
      )
      
      
      brep <- switch(
        estimator,
        "2SLS" = miive.2sls(d, g.b, r, est.only = TRUE, se, missing, var.cov),
        #"GMM"  = miive.gmm(d, g, r, est.only = TRUE), 
        # Not implemented
        # In other cases, raise an error
        stop(paste(
          "Invalid estimator:", 
          estimator, 
          "Valid estimators are: 2SLS")
        )
      )
      
      if (var.cov){
        
        v.b <- estVarCovarCoef(bsample, g.b, brep$eqn, pt, ordered)
        
        c(brep$coefficients, v.b$coefficients)
        
      } else {
        
        brep$coefficients
        
      }
      
    }, bootstrap)
    
    # Replace the estimated variances of estimates 
    # with the boostrap estimates
    boot.mat <- boot.results$t[
      stats::complete.cases(boot.results$t),
    ]
    
    results$bootstrap.true <- nrow(boot.mat)
    
    # Need to determine which columns are constants
    if (results$bootstrap.true > 0){
      
      nearZero <- apply(boot.mat , 2, stats::var, na.rm=TRUE) < 1e-16
      
    } else {
      
      nearZero <- FALSE
      
    }
    
    coefCov <- stats::cov(boot.mat)
    
    if(any(nearZero)) {
      coefCov[nearZero,nearZero] <- round(
        coefCov[nearZero,nearZero], 1e-10
      )
    }
    
    rownames(coefCov) <- colnames(coefCov) <- c(
      names(results$coefficients), 
      if (var.cov) { paste0(names(v$coefficients))} else { NULL }
    )
    
    results$eqn <- lapply(results$eqn, function(eq) {
      eq$coefCov <- coefCov[
          names(eq$coefficients),
          names(eq$coefficients),
          drop = FALSE
      ]
      eq
    })
    
    results$coefCov <- coefCov[
      names(results$coefficients),
      names(results$coefficients),
      drop = FALSE
    ]
    
    if (var.cov){
      v$coefCov <- coefCov[
        names(v$coefficients),
        names(v$coefficients),
        drop = FALSE
      ]
      
    }
    # Store the boot object as a part of the result object. 
    # This is useful for calculating CIs or
    # other bootstrap postestimation.
    results$boot <- boot.results

  } # end bootstrap
  
  
  # Disable twostep procedure in the current release
  
  if (missing == "twostage" & se == "standard" & var.cov == TRUE){

    twoStageCoefCov <- estTwoStageVarCovSE(g,v,results$eqn,pt)

    v$coefCov <- twoStageCoefCov[
       names(v$coefficients),
       names(v$coefficients),
       drop = FALSE
    ]
  }
  
  # assemble return object
  results$model          <- model
  results$estimator      <- estimator
  results$se             <- se
  results$var.cov        <- var.cov
  results$missing        <- missing
  results$boot.ci        <- boot.ci
  results$bootstrap      <- bootstrap
  results$call           <- match.call()
  results$ordered        <- ordered
  results$pt             <- unclass(pt)
  results$eqn.unid       <- d.un
  results$r              <- r
  results$v              <- v

  class(results)  <- "miive"
  
  results
}

