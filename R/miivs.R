#' Model-implied instrumental variable (MIIV) search 
#'
#' A key step in the MIIV-2SLS approach is to transform the SEM by 
#' replacing the latent variables with their scaling indicators 
#' minus their errors.  Upon substitution the SEM is transformed 
#' from a model with latent variables to one containing observed 
#' variables with composite errors. The miivs function 
#' automatically makes this transformation. The miivs function
#' will also identify equation-specific model-implied 
#' instrumental variables in simultaneous equation models
#' without latent variables. 
#'
#' @param model A model specified using lavaan model syntax. See the 
#' \code{model} argument within the \code{\link[lavaan]{lavaanify}} 
#' function for more information. See the documentation below for a 
#' description of how to specify the scaling  indicator in latent 
#' variable models and impose equality constraints on the parameter 
#' estimates. 
#'  
#' @details 
#' 
#' \itemize{
#' \item{\code{model}} {
#' 
#'   A model specified using the model syntax employed by \pkg{lavaan}. 
#'   The following model syntax operators are currently supported: =~,
#'   ~, ~~ and *. See below for details on default behaviors, 
#'   how to specify the scaling indicator in latent variable models, 
#'   and how to impose equality constraints on the parameter estimates. 
#'   
#'   \strong{Example using Syntax Operators}
#'   
#'   In the model below, 'L1 =~ Z1 + Z2 + Z3'  indicates the 
#'   latent variable L1 is measured by 3 indicators, Z1, Z2, and Z3. Likewise,
#'   L2 is measured by 3 indicators, Z4, Z5, and Z6. The statement
#'   'L1 ~ L2' specifies latent variable L1 is regressed on latent 
#'   variable L2. 'Z1 ~~ Z2' indicates the error of Z2 is allowed to 
#'   covary with the error of Z3. The label
#'   LA3 appended to Z3 and Z6 in the measurement model equations 
#'   constrains the factor loadings for Z3 and Z6 to equality. For 
#'   additional details on constraints see Equality Constraints 
#'   and Parameter Restrictions.
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
#'   Following the \pkg{lavaan} model syntax, latent variables are defined 
#'   using the \code{=~} operator.  For first order factors, the scaling 
#'   indicator chosen is the first observed variable on the RHS of an 
#'   equation. For the model  below \code{Z1} would be chosen as the 
#'   scaling indicator for \code{L1} and \code{Z4} would be chosen as 
#'   the scaling indicator for \code{L2}.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + Z2 + Z3
#'      L2 =~ Z4 + Z5 + Z6
#'   '}
#'   
#'   
#'   \strong{Equality Constraints and Parameter Restrictions}
#'   
#'   Within- and across-equation equality constraints on the factor loading
#'   and regression coefficients can be imposed directly in the model syntax. 
#'   To specify equality constraints between different parameters equivalent
#'   labels should be prepended to the variable name using the 
#'   * operator. For example, we could constrain the factor 
#'   loadings for two non-scaling indicators of latent factor \code{L1} to 
#'   equality using the following  model syntax.
#'   
#'   \preformatted{model <- '
#'      L1 =~ Z1 + LA2*Z2 + LA2*Z3
#'      L2 =~ Z4 + Z5 + Z6
#'   '}
#'   
#'   Researchers can also constrain the factor loadings and regression 
#'   coefficients to specific numeric values in a similar fashion. Below 
#'   we constrain the regression coefficient  of \code{L1} on \code{L2} 
#'   to \code{1}.
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
#'   indicator that would have been assigned to the first lower-order 
#'   factor \code{L1}. The intercepts for lower-order latent variables 
#'   are set to zero, by default
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
#'   \pkg{MIIVsem} will automatically include the intercepts of any 
#'   observed or latent endogenous variable. The intercepts
#'   for any scaling indicators and lower-order latent variables are
#'   set to zero. Covariances among exogenous latent
#'   and observed  variables are included by default when \code{
#'   var.cov = TRUE}. Where appropriate the covariances of the errors
#'   of latent and observed dependent variables are also automatically 
#'   included in the model specification. These defaults correspond 
#'   to those used by \pkg{lavaan} and \code{auto = TRUE}, except that 
#'   endogenous latent variable intercepts are estimated by default, 
#'   and the intercepts of scaling indicators are fixed to zero.
#'   
#'   \strong{Invalid Specifications}
#'   
#'   Certain model specifications are not currently supported.  For example,
#'   the scaling indicator of a latent variable is not permitted to
#'   cross-load on another latent variable. In the model below
#'   \code{Z1}, the scaling indicator for L1, cross-loads on the latent 
#'   variable \code{L2}. Executing a search on the model below will 
#'   result in the warning: \emph{miivs: scaling indicators with a factor 
#'   complexity greater than 1 are not currently supported}.
#'   
#'   \preformatted{model <- '
#'     L1 =~ Z1 + Z2 + Z3
#'     L2 =~ Z4 + Z5 + Z6 + Z1
#'   '}
#'   
#'   In addition, \pkg{MIIVsem} does not currently support relations
#'   where the scaling indicator of a latent variable is also the 
#'   dependent variable in a regression equation.  For example, the
#'   model below would not be valid under the current algorithm.
#'   
#'   \preformatted{model <- '
#'     L1 =~ Z1 + Z2 + Z3
#'     Z1  ~ Z4
#'     Z4  ~ Z5 + Z6
#'   '}
#'    
#'   
#'   }
#' }
#' 
#' The \code{miivs} function displays a table containing the following 
#' information for each equation in the system:
#' 
#' \itemize{
#'  \item \code{LHS} The "dependent" variable.
#'  \item \code{RHS} The right hand side variables of the transformed equation.
#'  \item \code{MIIVs} The model-implied instrumental variables for each equation.
#' }
#' 
#' @return A list of model equations.
#'
#' @references 
#' 
#' Bollen, K. A. (1996).	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 
#' Bentler, P. M., and Weeks, D. G. (1980). Linear Structural 
#' Equations with Latent Variables. \emph{Psychometrika}, 45, 289–308.                 
#' 	
#' @example example/bollen1989-miivs.R
#'
#' @seealso \code{\link{miive}}
#' 
#' @export
miivs <- function(model){

  pt <- lavaan::lavaanify( model, 
                           auto = TRUE, 
                           meanstructure = TRUE )
  

  #------------------------------------------------------------------------#
  # Parse parTable and add any equality or numeric constraints .           
  #
  #   More complicated constraints are possible during estimation but
  #   would require a more systematic method for parsing the parTable.  
  #------------------------------------------------------------------------#
  bad.ops <- c("==", "~*~","<~")
  
  # Stop if there are any constraints entered using the '==' operator.
  if (length(pt$lhs[pt$op %in% bad.ops & pt$user != 2]) > 0) {
    stop(paste("miivs: MIIVsem does not currently support",
               "the following operators:", 
               paste0(bad.ops,collapse = ", "),"."))
  }

  # Store any usable constraints in the mlabel column of pt. By using
  # the label column we take care of any labels entered using the "*" op.
  pt$mlabel <- pt$label
  
  # Numeric contraints eneterd using '*' op. in the model syntax should
  # also be added to the mlabel column
  condNum <- !(pt$op == "=~" & !duplicated(pt$lhs)) & 
             (!is.na(pt$ustart) & pt$free == 0)
  if (length(pt$ustart[condNum]) > 0){
    pt[condNum,]$mlabel <- pt[condNum,]$ustart
  }
  

  tmpMarkers <- pt[pt$op == "=~",]$rhs[
    which(!duplicated(pt[pt$op == "=~",]$lhs))
  ]
  
  # For now, throw an error if the scaling indicator cross loads
  if (length(tmpMarkers) > 0){
    for(i in 1:length(tmpMarkers)){
      if(length(pt[pt$op == "=~" & pt$rhs == tmpMarkers[i],"lhs"]) > 1){
        stop(paste("miivs: scaling indicators with a factor complexity", 
                   "greater than 1 are not currently supported."))
      }
      if(length(pt[pt$op == "~" & pt$lhs == tmpMarkers[i],"lhs"]) > 0){
        stop(paste("miivs: scaling indicators cannot be depdendent", 
                   "variables in regression equations."))
      }
    }
    pt[pt$op == "=~" & pt$rhs %in% tmpMarkers,]$mlabel <- NA
  }
  
  pt$mlabel[pt$mlabel == ""] <- NA
  
  #------------------------------------------------------------------------#
  # Build the Bentler-Weeks model matrices based on the parTable.          #
  #------------------------------------------------------------------------#
  
  latVars <- unique(pt$lhs[pt$op == "=~"])
  
  obsVars <- setdiff(
    unique(c(pt$lhs[pt$op != "=="],pt$rhs[pt$op != "=="])), 
    latVars
  )
  
  # Indicators of latent variables and DVs in regressions
  endVars <- unique(c(pt$rhs[pt$op == "=~"], pt$lhs[pt$op == "~"]))  
  
  errVars <- paste("e.",endVars,sep="")
  
  # Exogenous and error terms of endogenous variables
  exoVars <- c(
    setdiff(
      unique(c(pt$rhs[pt$op!="=="], pt$lhs[pt$op!="=="])), 
      endVars
    ), errVars)
  
  allVars <- c(endVars,exoVars)
  
  # Matrix orders
  n <- length(exoVars)
  m <- length(endVars)
  s <- m + n
  
  # Model matrices
  
  # paths between dependent and independent variables
  gamma <-  matrix(0, m, n, dimnames = list(endVars, exoVars)) 
  
  # paths between dependent variables
  beta <- matrix(0, m, m, dimnames = list(endVars, endVars)) 
  
  # covariances
  Phi <-   matrix(0, n, n, dimnames = list(exoVars, exoVars))   
  
  # Populate the matrices setting free parameters to NA and 
  # others to their fixed values. 
  # Edit:  here we only fix the scaling indicators for now.
  paramValues <- ifelse(
    pt$free == 0 & is.na(pt$mlabel),
    pt$ustart,
    NA
  )
  
  # Fill gamma matrix
  for(i in which(pt$op == "=~" & pt$lhs %in% exoVars)){
    
    gamma[pt$rhs[i],pt$lhs[i]] <- paramValues[i]
    
  }
  
  for(i in which(pt$op == "~" & pt$rhs %in% exoVars)){
    
    gamma[pt$lhs[i],pt$rhs[i]] <- paramValues[i]
  
  }
  
  # Error terms have coefficient 1 by definition
  gamma[,(n-m+1):n] <- diag(m)
  
  
  # Fill beta matrix
  for(i in which(pt$op == "=~" & ! pt$lhs %in% exoVars)){
    
    beta[pt$rhs[i],pt$lhs[i]] <- paramValues[i]
  }
  
  for(i in which(pt$op == "~" & ! pt$rhs %in% exoVars)){
    
    beta[pt$lhs[i],pt$rhs[i]] <- paramValues[i]
    
  }
  
  # Covariances between endogenous LVs are actually 
  # covariances between error terms
  lhs <- ifelse(
    pt$lhs %in% endVars,
    paste("e.",pt$lhs,sep=""),
    pt$lhs
  
    )
  rhs <- ifelse(
    pt$rhs %in% endVars,
    paste("e.",pt$rhs,sep=""),
    pt$rhs
  )
  
  for(i in which(pt$op == "~~")){
    
    Phi[lhs[i],rhs[i]] <- Phi[rhs[i],lhs[i]] <- paramValues[i]
    
  }
  
  #------------------------------------------------------------------------#
  # Calculate full Sigma (without selection matrix G)                      #
  #------------------------------------------------------------------------#
  
  Beta  <- matrix(0,s,s, dimnames = list(allVars, allVars))
  Gamma <- matrix(0,s,n, dimnames = list(allVars, exoVars))
  I     <- diag(nrow(Beta))
  
  # Populate the matrices
  Gamma[1:m,] <- gamma
  Gamma[(m+1):s,] <- diag(n)
  Beta[1:m,1:m]<-beta
  
  # Calculate the model implied covariance matrix. We need to define 
  # a custom matrix operator. Normally, anything multiplied by NA 
  # (free parameters in our case) is NA. However we need 0*NA = 0. 
  `%naproduct%` <- function(x, y) {
    as.matrix(
      Matrix::Matrix(x, sparse = T) %*% 
      Matrix::Matrix(y, sparse = T)
    )
  }
  
  # Sigma depends on solve(I-Beta), which gives total effects between 
  # endogenous variables. However, this cannot be calculated when 
  # NAs are present, so we need to do it "by hand"
  
  BetaI <- I-Beta
  BetaI[is.na(BetaI)] <- 0 # First set NAs to zero to be able to solve
  BetaI <- solve(BetaI)
  
  # Then solve setting all non-zero effects to ones to find out element 
  # that are influenced by NAs.
  BetaNA <- I-(is.na(Beta) | Beta != 0)
  
  # Try to invert BetaNA
  
  trySolve <- function(mat){
    class(try(solve(mat),silent=T))=="matrix"
  }
  if (trySolve(BetaNA)){ 
    BetaNA <- solve(BetaNA) 
  } else {
    # fill the matrix with noise
    nz <- length(BetaNA[BetaNA==-1])
    BetaNA[BetaNA==-1] <- stats::runif(nz)
    BetaNA <- solve(BetaNA) 
  }
  BetaI[BetaNA != 0 & BetaI == 0] <- NA
  
  Sigma <- BetaI %naproduct% Gamma %naproduct% 
           Phi %naproduct% t(Gamma) %naproduct% t(BetaI)
  
  #------------------------------------------------------------------------#
  # Identify MIIVs for each equation in the system.                        #
  #------------------------------------------------------------------------#
  
  # Matrix of all regressions
  gamBeta <- cbind(gamma,beta)
  
  # Build the MIIV models:Loop over all variables that receive a path 
  eqns <- list()
  
  for(dv in unique(rownames(gamBeta)[which(apply(is.na(gamBeta),1,any))])){
    
    eq <- list()
    
    # Get anything on the rhs of dv, including composite error terms
    vars <- c(dv,colnames(gamBeta)[c(which(gamBeta[dv,]!=0), which(is.na(gamBeta[dv,])))])
    
    # Add DVs and IVs before potential scaling indicator substitution. 
    eq$EQnum <- NA
    eq$EQmod <- NA
    eq$DVobs <- NA
    eq$IVobs <- NA
    eq$DVlat <- dv
    eq$IVlat <- setdiff(vars[-1], errVars)

    compositeDisturbance <- paste("e.",vars[1],sep="")
    
    # 1) Choose marker indicators for each latent variable
    
    markers  <- NULL
    
    # For each variable on the RHS including error terms
    for(var in vars){
      
      if(var %in% latVars){
        
        # Choose the first indicator with loading fixed to 1
        marker <- rownames(gamBeta)[which(gamBeta[,var]==1)[1]]
        
        compositeDisturbance <- c(
          compositeDisturbance,paste("e.",marker,sep="")
        )
        
        # Deal with higher order factors
        while(marker %in% latVars){
          
          eq$EQmod <- "measurement"
          
          marker <- rownames(gamBeta)[which(gamBeta[,marker]==1)[1]]
          
          compositeDisturbance <- c(
            compositeDisturbance,
            paste("e.",marker,sep="")
          )
          
        }
        
        markers <- c(markers, marker)
        vars[vars == var] <- marker
        
      }
      
    }
    
    eq$DVobs   <- vars[1]
    eq$IVobs   <- setdiff(vars[-1], errVars)
    eq$CDist   <- compositeDisturbance
    eq$MIIVs   <- NA
    eq$markers <- eq$markers
    eqns <- c(eqns,list(eq))
    
  }
  
  # FIXME: 
  ## are there any scaling indicators listed as separate DVs
  # dups  <- duplicated(lapply(eqns,"[[","DVobs"))
  # 
  # if (any(dups)){
  #   
  #   cond   <- unlist(lapply(eqns,function(x){x$DVobs %in% x$IVobs}))
  #   eqns.x <- eqns[ cond] # regression relationships
  #   eqns   <- eqns[!cond] # measurement relationships
  # 
  #   
  #   dup.dvs <- unlist(lapply(eqns.x,"[[","DVobs"))
  #   
  #   eqns <- lapply(eqns, function(eq) {
  #     
  #     idx <- which(dup.dvs %in% eq$DVobs)
  #     
  #     if (length(idx) > 0){
  #       
  #       if (length(intersect(eq$IVlat, eqns.x[[idx]]$IVlat))>0){
  #         stop("miivs: miivs cannot parse this model syntax.")
  #       }
  #       
  #       non.dup.obs <- eqns.x[[idx]]$IVobs[
  #         !eqns.x[[idx]]$IVobs %in% eqns.x[[idx]]$DVobs  
  #       ]
  #       
  #       non.dup.lat <- eqns.x[[idx]]$IVlat[
  #         !eqns.x[[idx]]$IVobs %in% eqns.x[[idx]]$DVobs 
  #       ]
  #       
  #       eq$EQmod <- "mixed"
  #       eq$VRtyp <- c(rep(T,length(eq$IVlat)), rep(F,length(non.dup.obs)))
  #       eq$IVlat <- c(eq$IVlat, non.dup.lat)
  #       eq$IVobs <- c(eq$IVobs, non.dup.obs)
  #       eq$CDist <- unique(c(eq$CDist, eqns.x[[idx]]$CDist))
  # 
  #     }
  #     
  #     eq
  #     
  #   })
  #   
  # }
  
    
  for (j in 1:length(eqns)){
    
    eqns[[j]]$EQnum <- j
    # 2) Take a submatrix Sigma_e, choosing n columns of S that correspond 
    #    to the error term of the dependent variable and error terms of the 
    #    marker indicators and rows that correspond to k observed variables 

    Sigma_e <- Sigma[,eqns[[j]]$CDist, drop = FALSE]
    
    # 3) Collapse matrix Sigma_e as k vector where an element receives 1 if all 
    #    elements of  Sigma_e are zero on a corresponding row. This vector 
    #    indicates which observed variables are uncorrelated with all 
    #    components fo the composite error term (vector e).
    
    e <- apply(Sigma_e == 0 & ! is.na(Sigma_e),1,all)
    
    # 4) Take a submatrix Sigma_i, choosing (n-1) columns of S that correspond 
    #    to the marker indicators and rows that correspond to k observed 
    #    variables.
    
    Sigma_i <- Sigma[,eqns[[j]]$markers, drop = FALSE]
    
    # 5) Collapse matrix Sigma_i as k vector where an element receives 1 if 
    #    any elements of  Sigma_e is nonzero on a corresponding row. This 
    #    vector indicates which observed variables have at least one non-zero 
    #    implied correlation with any of the marker indicators (vector i).
    
    i <- apply(Sigma_i != 0 | is.na(Sigma_i),1,all)
    
    # 6) Do an element wise multiplication of i and e. This produces a 
    #    vector indicating which indicators have zero implied correlations 
    #    between the composite error term and at least one implied non-zero 
    #    correlation with the marker indicators and thus contains the indices 
    #    of all valid instruments.
    
    eqns[[j]]$MIIVs <- names(which((e&i)[obsVars]))
    
    # Add equation type: we did this for higher order factors above. 
    eqns[[j]]$EQmod <- if(is.na(eqns[[j]]$EQmod)){
      if (eqns[[j]]$DVlat %in% pt$rhs[pt$op =="=~"] ){
        "measurement"
      } else {
        "regression"
      }
    } else {
      eqns[[j]]$EQmod 
    }
    
  }
  
  eqns <- lapply(eqns,function(eq){eq$markers <- NULL; eq})
  
  #------------------------------------------------------------------------#
  # Add the labels
  #    Better to do this at the end in case of constraints among
  #    higher order factors. 
  #------------------------------------------------------------------------#
  
  lhsi <- c(pt$rhs[pt$op == "=~" & duplicated(pt$lhs)],
            pt$lhs[pt$op == "~"])
  
  rhsi <- c(pt$lhs[pt$op == "=~" & duplicated(pt$lhs)],
            pt$rhs[pt$op == "~"])
  
  labi <- c(pt$mlabel[pt$op == "=~" & duplicated(pt$lhs)],
            pt$mlabel[pt$op == "~"])
  
  for ( j in 1:length(eqns)){ 
    td <- eqns[[j]]$DVlat
    ti <- eqns[[j]]$IVlat
    eqns[[j]]$Label <- labi[which(lhsi %in% td)][
      match(ti, rhsi[which(lhsi %in% td)])
    ]
  }
  
  res <- list(
    eqns = eqns, 
    pt = pt,
    matrices = 
      list(
        Sigma = Sigma, 
        Beta = Beta, 
        BetaI = BetaI, 
        Gamma = Gamma, 
        Phi = Phi
      )
  )
  
  class(res) <- "miivs"
  res
  
}
