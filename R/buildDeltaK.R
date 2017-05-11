#'@keywords internal

  buildDeltaK <- function(eqns, s.cov, s.mean,acov.names){
    #==============================================================#
    K <- do.call("rbind",
      lapply(eqns, function(eq){
        Svv <- s.cov[eq$MIIVs, eq$MIIVs, drop = FALSE]
        Svz <- s.cov[eq$MIIVs, eq$IVobs, drop = FALSE]
        Svy <- s.cov[eq$MIIVs, eq$DVobs, drop = FALSE]
        muV <- s.mean[eq$MIIVs]
        mu  <- c(s.mean[eq$DVobs], s.mean[eq$IVobs])
        B   <- eq$coefficients[-1]
        
        derivMu     <- function(mu, B){
          mu[1] - t(mu[-1])%*%B
        }
        dMu.numeric <- t(numDeriv::jacobian(derivMu,mu,B=B))
        dMu.names <- cbind(names(mu), "1")
        
        #==============================================================#
        derivSvv     <- function(vechSvv, Svz, Svy){
          Svv <- revVech(vechSvv)
          solve(t(Svz)%*%solve(Svv)%*%Svz) %*% t(Svz)%*%solve(Svv) %*% Svy
        }
        
        dSvv.numeric <- t(numDeriv::jacobian(derivSvv, 
                                             Svv[lower.tri(Svv, diag=TRUE)], 
                                             Svz = Svz,
                                             Svy = Svy))
        names.cov <- outer(colnames(Svv), colnames(Svv), function(x, y){
          paste0(x, "~~",y)
        })
        names.cov <- names.cov[lower.tri(names.cov, diag = TRUE)]
        dSvv.names <- do.call("rbind",strsplit(names.cov,"~~"))
        #==============================================================#
        derivSvz     <- function(Svz, Svv, Svy){
          U  <- solve(t(Svz) %*% solve(Svv) %*% Svz)
          U %*% (t(Svz)%*% solve(Svv) %*% Svy)
        }
        dSvz.numeric <- t(numDeriv::jacobian(derivSvz, 
                                             Svz, 
                                             Svv = Svv, 
                                             Svy = Svy))
        dSvz.names <- as.matrix(expand.grid(rownames(Svz),colnames(Svz)))
        #==============================================================#
        derivSvy     <- function(Svy, Svz, Svv){
          solve(t(Svz)%*%solve(Svv)%*%Svz) %*% t(Svz)%*%solve(Svv) %*% Svy
        }
        dSvy.numeric <- t(numDeriv::jacobian(derivSvy, 
                                             Svy, 
                                             Svz = Svz, 
                                             Svv = Svv))
        dSvy.names <- as.matrix(expand.grid(rownames(Svy),colnames(Svy)))
  
        #==============================================================#
        covDeriv <- rbind(dSvv.numeric, dSvz.numeric, dSvy.numeric)
        covNames <- rbind(dSvv.names, dSvz.names, dSvy.names)
        muDeriv  <- rbind(dMu.numeric)
        muNames  <- rbind(dMu.names)
        covVec <- rep(NA, nrow(covNames ))
        
        for(j in 1:nrow(covNames)){ 
          
          covVec[j] <- match(
            paste0(covNames[j,1],"~~",covNames[j,2]),
            acov.names
          )
          
          if (is.na(covVec[j])) {
            covVec[j] <- match(
              paste0(covNames[j,2],"~~", covNames[j,1]),
              acov.names
            )
          }
          
        }
        
        muVec <- rep(NA, nrow(muNames))
        muVec <- match(paste0(muNames[,1],"~",muNames[,2]),acov.names)
        
        row.k             <- matrix(0, ncol(covDeriv)+1, length(acov.names))
        row.k[ 1, muVec ] <- t(muDeriv)
        row.k[-1, covVec] <- t(covDeriv)
        row.k
      })
    )
    #==============================================================#
    return(K)
  }
