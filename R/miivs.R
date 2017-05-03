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
#' documentation below for a description of how to specify the scaling 
#' indicator in latent variable models and impose equality constraints on the
#' parameter estimates. 
#' 
#' @section Scaling Indicators:
#' Following the \code{lavaan} model syntax, latent variables are defined 
#' using the \code{"=~"} operator.  For first order factors, the scaling indicator 
#' chosen is the first observed variable on the RHS of an equation. For the model 
#' below \code{Z1} would be chosen as the scaling indicator for \code{L1} and  
#' \code{Z4} would be chosen as the scaling indicator for \code{L2}.
#' 
#' \preformatted{model <- '
#'     L1 =~ Z1 + Z2 + Z3
#'     L2 =~ Z4 + Z5 + Z6
#' '
#' }
#' 
#' @section Higher-order Factor Models:
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
#' @return A list of model equations.
#'
#' @references 
#' 
#' Bollen, K. A. 1996.	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 	
#' @example example/bollen1989-miivs.R
#'
#' @seealso \code{\link{miive}}
#' 
#' @export

miivs <- function(model){

  pt <- lavaan::lavaanify(model, auto = TRUE)
  
  # Remove any covariance or variance constrained to zero
  var.covar.zero <- which(
    pt$op     == "~~" & 
    pt$ustart == 0
  )
  
  if (length(var.covar.zero) > 1){
    pt <- pt[-var.covar.zero, , drop = FALSE]
  }

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
  # also be added to the mlabel column.
  condNum <- !(pt$op == "=~" & !duplicated(pt$lhs)) & !is.na(pt$ustart)
  if (length(pt$ustart[condNum]) > 0){
    pt[condNum,]$mlabel <- pt[condNum,]$ustart
  }
  
  # Remove any scaling indicators
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
    }
    pt[pt$op == "=~" & pt$rhs %in% tmpMarkers,]$mlabel <- NA
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
  
  # Exogenous and error terms of endogenous variables
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
    as.matrix(Matrix::Matrix(x, sparse = T) %*% Matrix::Matrix(y, sparse = T))
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
    BetaNA[BetaNA==-1] <- runif(nz)
    BetaNA <- solve(BetaNA) 
  }
  BetaI[BetaNA != 0 & BetaI == 0] <- NA
  
  Sigma <- BetaI %naproduct% Gamma %naproduct% Phi %naproduct% 
           t(Gamma) %naproduct% t(BetaI)
  
  #------------------------------------------------------------------------#
  # Identify MIIVs for each equation in the system.                        #
  #------------------------------------------------------------------------#
  
  # Matrix of all regressions
  gamBeta <- cbind(gamma,beta)
  
  # Build the MIIV models:Loop over all variables that receive a path 
  eqns <- list(); cntr <- 1;
  
  for(dv in unique(rownames(gamBeta)[which(apply(is.na(gamBeta),1,any))])){
    
    eq <- list()
    
    # Get anything on the rhs of dv, including composite error terms
    vars <- c(dv,colnames(gamBeta)[c(which(gamBeta[dv,]!=0), which(is.na(gamBeta[dv,])))])
    
    # Add DVs and IVs before potential scaling indicator substitution. 
    eq$EQnum <- cntr
    eq$EQmod <- NA
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
        compositeDisturbance <- c(compositeDisturbance,paste("e.",marker,sep=""))
        
        # Deal with higher order factors
        while(marker %in% latVars){
          
          eq$EQmod <- "measurement"
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
    
    # Add equation type: we did this for higher order factors above. 
    eq$EQmod <- if(is.na(eq$EQmod)){
      if ( (eq$IVlat != eq$IVobs) && (eq$DVlat == eq$DVobs) ){
        "measurement"
      } else {
        "structural"
      }
    }
    
    cntr <- cntr + 1
    
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
  
  res <- list(
    eqns = eqns, 
    pt = pt, 
    matrices = list(
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
