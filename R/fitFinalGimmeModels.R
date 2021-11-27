#'@keywords internal
fitFinalGimmeModels <- function(ts_list_obs, meas_model, lv_model, miiv.dir, lv_final_estimator){
  
  
   # ts_list_obs = dat$lvgimme$ts_list_obs
   # meas_model  = dat$lvgimme$model_list_dfa
   # lv_model    = lapply(store$syntax, function(x){x[!grepl("0\\*", x)]})
   # miiv.dir    = file.path(dat$out,"miiv")
   # lv_final_estimator = lv_final_estimator

  ts_list_obs <- lapply(ts_list_obs, function(df){
    
    first           <- df[1:(nrow(df)-1), ] 
    second          <- df[2:(nrow(df)  ), ]
    ts_lc           <- data.frame(first, second)
    colnames(ts_lc) <- c(paste0(colnames(df), "lag"), colnames(df))
    ts_lc
      
  })
  
  
  df_all <- do.call("rbind",lapply(seq_along(ts_list_obs), function(i){
    
    pt  <- lavaan::lavParTable(c(meas_model[[i]], lv_model[[i]]),fixed.x=FALSE)
    
    # remove any nonsense paths that have been fixed to zero
    pt <- pt[pt$free != 0,      ]
    pt <- pt[pt$op   != ".==.", ]
    pt <- pt[!(pt$op   == "~~" & pt$lhs == pt$rhs), ]
    pt <- pt[pt$op != "~1",]
    pt <- pt[!(grepl("lag", pt$lhs) & pt$op == "~~"),]
    
    mod <- unlist(lapply(seq_along(1:nrow(pt)), function(j){
      pt <- pt[j,]
      if(pt$label == ""){
        paste0(pt$lhs, pt$op, pt$rhs)
      } else {
        paste0(pt$lhs, pt$op,pt$label,"*", pt$rhs)
      }
    }))
    
    
    if(lv_final_estimator == "miiv"){
      
      so <- miivs(mod)
      
      so$eqns <- lapply(so$eqns, function(x){
        
        if(x$EQmod == "measurement"){
          
          #lv <- ifelse(grepl("lag", x$IVlat), x$IVlat, paste0(x$IVlat,"lag"))
          lv <- x$IVlat
          non_scaling_indicators <- pt[pt$op=="=~" & pt$lhs %in% lv & pt$label !="", "rhs"]
          x$MIIVs <- setdiff(non_scaling_indicators, c(x$DVobs, x$IVobs))
          
        } else if (x$EQmod == "regression"){
          
          rhs_lv  <- x$IVlat
          x$MIIVs <- pt[pt$op=="=~" & pt$lhs %in% rhs_lv & pt$label !="", "rhs"]
          
        }
        x
      })
      
      df <- estimatesTable(miive(so, ts_list_obs[[i]], overid.degree = 1, overid.method = "stepwise.R2"), sarg = TRUE)
      df <- df[df$op != "~1",]
      
      
    } else if (lv_final_estimator == "pml"){
      
      df <- lavaan::parameterEstimates(lavaan::sem(mod, ts_list_obs[[i]]))
      df <- df[df$op != "~1",]
            
    }
    
    subj <- names(ts_list_obs)[i]
    
    cbind(subj ,df)
  }))
  

  if(!is.null(miiv.dir)){
    dir.create(miiv.dir,showWarnings = FALSE)
    
    utils::write.csv(df_all,file.path(miiv.dir, "indPathEstimates.csv") ,row.names = FALSE)
  }
  
  return(df_all)
  
}

