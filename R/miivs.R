#' Model-implied instrumental variable (MIIV) search 
#'
#' A key step in the MIIV-2SLS approach is to transform the SEM by replacing the 
#' latent variables with their scaling indicators minus their errors.  Upon 
#' substitution the SEM is transformed from a model with latent variables to 
#' one containing observed variables with composite errors.  The miivs function 
#' automatically makes this transformation.
#'
#' @param model A model specified using lavaan model syntax. See the \code{model} 
#' argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' A model specified using lavaan model syntax. See the \code{model} argument within 
#' the \code{\link[lavaan]{lavaanify}} function for more information. See the
#' documentation below for a description of how to specifying the scaling 
#' indicator in latent variable models and impose equality constraints on the
#' parameters. 
#' 
#' @param miivs.out A logical indicating whether or not to print the MIIVs as an 
#' object for later use. Default is \code{FALSE}.
#'
#' 
#' @section Scaling Indicators:
#' Following the lavaan model syntax, latent variables are defined using the
#' \code{=~} operator.  For first order factors, the scaling indicator chosen
#' is the first observed variable on the RHS of an equation. For the model below
#' \code{Z1} would be chosen as the scaling indicator for \code{L1} and  
#' \code{Z4} would be chosen as the scaling indicator for \code{L2}.
#' 
#' \preformatted{model <- '
#'     L1 =~ Z1 + Z2 + Z3
#'     L2 =~ Z4 + Z5 + Z6
#' '
#' }
#' 
#' For higher-order factor models, the scaling indicator for the higher order
#' factor is taken from the scaling indicator that would have been 
#' assigned to the first lower-order factor on the RHS of the equation. 
#' For example, in the model below, the  scaling indicator for the higher-order 
#' factor \code{H1} is taken to be \code{Z1}, the scaling indicator 
#' that would have been assigned to the first lower-order factor \code{L1}.
#' 
#'
#' \preformatted{model <- '
#'     H1 =~ L1 + L2 
#'     L1 =~ Z1 + Z2 + Z3
#'     L2 =~ Z4 + Z5 + Z6
#' '
#' }
#' 
#' @section Equality Constraints:
#' 
#' Within- and across-equation equality constraints on the factor loading
#' and regression coefficients can be imposed directly in the model syntax. 
#' To specify equality constraints between different parameters equivalent
#' labels should be prepended to the variable name using the \code{*} operator.
#' For example, we could constrain the factor loadings for two non-scaling 
#' indicators of latent factor \code{L1} to equality using the following 
#' model syntax.  
#' 
#' \preformatted{model <- '
#'     L1 =~ Z1 + B1*Z2 + B1*Z3
#'     L2 =~ Z4 + Z5 + Z6
#' '
#' }
#' 
#' The factor loading and regression coefficients can also be constrained
#' to specific numberic values in a similar fashion. Below we constrain
#' the second factor loading for the \code{L1} factor to a value of \code{2}
#' 
#' \preformatted{model <- '
#'     L1 =~ Z1 + 2*Z2 + Z3
#'     L2 =~ Z4 + Z5 + Z6
#' '
#' }
#' 
#' @details 
#' 
#' The \code{miivs} function displays a table containing the following 
#' information for each equation in the system:
#' 
#' \itemize{
#'  \item \code{LHS} The "dependent" variable.
#'  \item \code{RHS} The right hand side variables of the transformed equation.
#'  \item \code{Composite Disturbance}  Elements of the composite errors in the transformed equation.
#'  \item \code{MIIVs} The model-implied instrumental variables for each equation.
#' }
#' 
#' @references 
#' 
#'  Bollen,	K. A. and	D. J.	Bauer.	2004.	Automating	the	Selection	of 
#' 	Model-Implied	Instrumental	Variables.	\emph{Sociological	Methods	and	
#' 	Research}, 32, 425-52.
#' 	
#' 	Bauldry, S.	2014.	miivfind: A command for identifying model-implied instrumental 
#' 	variables for structural equation models in Stata.	\emph{Stata Journal}, 14:4.
#' 	
#' 	Bentler, P. M., and Weeks, D. G. (1980). Linear structural equations with 
#' 	latent variables. Psychometrika, 45(3), 289–308.                 
#' 	
#' @example example/bollen1989a.R
#'   
#' @export

miivs <- function(model){

  
  pt <- lavaanify(model, auto = TRUE)
  
  #------------------------------------------------------------------------#
  # Parse parTable and add any equality or numeric constraints .           
  #
  #   More complicated constraints are possible during estimation but
  #   would require a more systematic method for parsing the parTable.  
  #------------------------------------------------------------------------#
  
  # Stop if there are any constraints entered using the '==' operator.
  if (length(pt$lhs[pt$op == "==" & pt$user != 2]) > 0) {
    stop(paste("miivs: MIIVsem does not currently support the '==' operator."))
  }

  # Store any usable constraints in the mlabel column of pt. By using
  # the label column we take care of any labels entered using the "*" op.
  pt$mlabel <- pt$label
  
  # Numeric contraints eneterd using '*' op. in the model syntax are 
  # also be added to the mlabel column.
  condNum <- !(pt$op == "=~" & !duplicated(pt$lhs)) & !is.na(pt$ustart)
  if (length(pt$ustart[condNum]) > 0){
    pt[condNum,]$mlabel <- pt[condNum,]$ustart
  }
  
  pt$mlabel[pt$mlabel == ""] <- NA
  
  #------------------------------------------------------------------------#
  # Build the Bentler-Weeks model matrices based on the parTable.          #
  #------------------------------------------------------------------------#
  
  latVars <- unique(pt$lhs[pt$op == "=~"])
  
  obsVars <- setdiff(unique(c(pt$lhs[pt$op != "=="], pt$rhs[pt$op != "=="])), latVars)
  
  # Indicators of latent variables and DVs in regressions
  endVars <- unique(c(pt$rhs[pt$op == "=~"], pt$lhs[pt$op == "~"]))  
  
  errVars <- paste("e.",endVars,sep="")
  
  # All variables, excluding exogenous and error terms of endogenous variables
  exoVars <- c(setdiff(unique(c(pt$rhs[pt$op!="=="], pt$lhs[pt$op!="=="])), endVars), errVars)
  
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
  paramValues <- ifelse(pt$free == 0 & is.na(pt$mlabel),pt$ustart,NA)
  
  
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
  lhs <- ifelse(pt$lhs %in% endVars,paste("e.",pt$lhs,sep=""),pt$lhs)
  rhs <- ifelse(pt$rhs %in% endVars,paste("e.",pt$rhs,sep=""),pt$rhs)
  
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
    r <- matrix(0,nrow(x),ncol(y))
    for(i in 1:nrow(r)){
      for(j in 1:ncol(r)){
        e <- x[i,]*y[,j]
        e[which(x[i,]==0)] <- 0
        e[which(y[,j]==0)] <- 0
        r[i,j] <- sum(e)
      }
    }
    rownames(r) <- rownames(x)
    colnames(r) <- colnames(y)
    return(r)
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
  BetaNA <- solve(BetaNA)
  BetaI[BetaNA != 0 & BetaI == 0] <- NA
  
  Sigma <- BetaI %naproduct% Gamma %naproduct% Phi %naproduct% 
           t(Gamma) %naproduct% t(BetaI)
  
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
    eq$DVlat <- dv
    eq$IVlat <- setdiff(vars[-1], errVars)

    compositeDisturbance <- paste("e.",vars[1],sep="")
    
    # 1) Choose marker indicators for each latent variable
    
    markers  <- NULL
    
    for(var in vars){
      
      if(var %in% latVars){
        
        # Choose the first indicator with loading fixed to 1
        marker <- rownames(gamBeta)[which(gamBeta[,var]==1)[1]]
        compositeDisturbance <- c(compositeDisturbance,paste("e.",marker,sep=""))
        
        # Deal with higher order factors
        while(marker %in% latVars){
          
          marker <- rownames(gamBeta)[which(gamBeta[,marker]==1)[1]]
          compositeDisturbance <- c(compositeDisturbance,paste("e.",marker,sep=""))
          
        }
        
        markers <- c(markers, marker)
        vars[vars == var] <- marker
        
      }
      
    }
  
    eq$DVobs <- vars[1]
    eq$IVobs <- setdiff(vars[-1], errVars)
    eq$CDist <- compositeDisturbance
    
    # 2) Take a submatrix Sigma_e, choosing n columns of S that correspond 
    #    to the error term of the dependent variable and error terms of the 
    #    marker indicators and rows that correspond to k observed variables 

    Sigma_e <- Sigma[,compositeDisturbance, drop = FALSE]
    
    # 3) Collapse matrix Sigma_e as k vector where an element receives 1 if all 
    #    elements of  Sigma_e are zero on a corresponding row. This vector 
    #    indicates which observed variables are uncorrelated with all 
    #    components fo the composite error term (vector e).
    
    e <- apply(Sigma_e == 0 & ! is.na(Sigma_e),1,all)
    
    # 4) Take a submatrix Sigma_i, choosing (n-1) columns of S that correspond 
    #    to the marker indicators and rows that correspond to k observed 
    #    variables.
    
    Sigma_i <- Sigma[,markers, drop = FALSE]
    
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
    
    eq$MIIVs <- names(which((e&i)[obsVars]))
    
    eqns <- c(eqns,list(eq))
  }
  
  #------------------------------------------------------------------------#
  # Add the labels
  #    Better to do this at the end in case of constraints among
  #    higher order factors. 
  #------------------------------------------------------------------------#
  lhsi <- c(pt$rhs[pt$op == "=~" & duplicated(pt$lhs)],pt$lhs[pt$op == "~"])
  rhsi <- c(pt$lhs[pt$op == "=~" & duplicated(pt$lhs)],pt$rhs[pt$op == "~"])
  labi <- c(pt$mlabel[pt$op == "=~" & duplicated(pt$lhs)],pt$mlabel[pt$op == "~"])
  
  for ( j in 1:length(eqns)){ # j <- 2
    td <- eqns[[j]]$DVlat
    ti <- eqns[[j]]$IVlat
    eqns[[j]]$Label <- labi[which(lhsi %in% td)][match(ti, rhsi[which(lhsi %in% td)])]
  }
  
  
  for (i in 1:length(eqns)){
    LHS <- paste(eqns[[i]]$DVobs, collapse = ", ")
    RHS <- paste(eqns[[i]]$IVobs, collapse = ", ")
    Instruments <- paste(eqns[[i]]$MIIVs, collapse = ", ")
    Disturbance <- paste("e.",eqns[[i]]$CDist, collapse = ", ", sep="")
    modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
    colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
    if (i == 1) {modeqns <- modtemp }
    if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
  }
  
  search <- list(eqns = eqns, df = modeqns, miivs.out = miivs.out)
  class(search) <- "miivs"
  search
}
