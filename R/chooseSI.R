#' identify best scaling indicator
#'@keywords internal
chooseSI <- function(model, data, scaling, lag = FALSE){
  
  if ( "miivs" == class(model) ){ 
    
    d  <- model$eqns 
    pt <- model$pt
    
  } else { 
    
    res <- miivs(model)
    d   <- res$eqns
    pt  <- res$pt
    
  } 
  
  latVars <- unique(pt$lhs[pt$op == "=~"])
  
  obsVars <- setdiff(
    unique(c(pt$lhs[pt$op != "=="],pt$rhs[pt$op != "=="])), 
    latVars
  )
  
  r2.list <- list(lv = "", r2 = "", cross.load = "", model = "")
  r2.list <- replicate(length(latVars), r2.list, simplify = FALSE)
  
  
  r2.list <- lapply(seq_along(r2.list), function(i){
    r2.list[[i]]$lv <- latVars[i]
    r2.list[[i]]$r2 <- getIndicatorR2s(pt[pt$op == "=~" & pt$lhs == latVars[i],"rhs"], data)
    r2.list[[i]]
  })
  
  
  indicator.list <- lapply(r2.list, function(x){names(x$r2)})
  
  r2.list <- lapply(seq_along(r2.list), function(i){
    
    common.indicators <- intersect(
        names(r2.list[[i]]$r2),
        unlist(indicator.list[-i])
      )
    
    if(length(common.indicators) > 0){
      r2.list[[i]]$cross.load <- common.indicators
    } else {
      r2.list[[i]]$cross.load <- NA
    }
    
    if(lag) {
      
      # inds <- names(sort(r2.list[[i]]$r2, decreasing = TRUE))
      # 
      # r2.list[[i]]$model <- paste0(
      #   r2.list[[i]]$lv, "=~", 
      #   paste0(names(sort(r2.list[[i]]$r2, decreasing = TRUE)), collapse ="+")
      # )
      # 
      # lag0 <- paste0(paste0(r2.list[[i]]$lv,"_L0") , "=~", paste0(inds,"*",paste0(inds,"_L0"), collapse ="+"))
      # 
      # lag1 <- paste0(paste0(r2.list[[i]]$lv,"_L1") , "=~", paste0(inds,"*",paste0(inds,"_L1"), collapse ="+"))
      # 
      # r2.list[[i]]$model <- c(lag0,lag1)

    } else {
      
      if(scaling != "first.indicator"){
        
        r2.list[[i]]$model <- paste0(
          r2.list[[i]]$lv, "=~", 
          paste0(names(sort(r2.list[[i]]$r2, decreasing = TRUE)), collapse ="+")
        )
      } else {
        r2.list[[i]]$model <- paste0(
          r2.list[[i]]$lv, "=~", 
          paste0(names(r2.list[[i]]$r2), collapse ="+")
        )
        
      }
      
    }


    r2.list[[i]]

  })

  return(r2.list)

}