#' provides model lists for latent variable gimme
#'@keywords internal
getModelList <- function(data, model, scaling = "first.indicator", lag = FALSE, dfm = FALSE){
  
  ind.list <- lapply(seq_along(data), function(i){
    lapply(seq_along(model[[i]]), function(j) { 
      chooseSI(model[[i]][[j]], data[[i]], scaling)
    })
  })
  
  if (scaling == "individual"){
    res <- lapply(seq_along(ind.list), function(i) {
      unlist(lapply(seq_along(ind.list[[i]]), function(j) { 
        inds <- paste0(names(sort(ind.list[[i]][[j]][[1]]$r2, decreasing = TRUE)))
        #inds_label <- paste0(inds,"*")
        inds_label <- paste0(ind.list[[i]][[j]][[1]]$lv,"_",inds,"*")
        inds_label[1] <- ""
        if(lag){
          lag0 <- paste0(paste0(ind.list[[i]][[j]][[1]]$lv) , "=~", paste0(inds_label,paste0(inds), collapse ="+"))
          lag1 <- paste0(paste0(ind.list[[i]][[j]][[1]]$lv,"lag") , "=~", paste0(inds_label,paste0(inds,"lag"), collapse ="+"))
          paste0(lag0,"\n",lag1, collapse = "")
        } else {
          paste0(paste0(ind.list[[i]][[j]][[1]]$lv) , "=~", paste0(inds_label,paste0(inds), collapse ="+"))
        }
        
      }))
    })
    
  } else if (scaling == "group"){
    
    res <- replicate(length(data), unlist(lapply( seq_along(ind.list[[1]]), function(i) {
      
      inds <- names(sort(colMeans(do.call("rbind", lapply( seq_along(ind.list), function(j) {
        ind.list[[j]][[i]][[1]]$r2
      }))), decreasing = TRUE))
      
      inds_label <- paste0(ind.list[[1]][[i]][[1]]$lv, "_", inds,"*")
      inds_label[1] <- ""
      
      if(lag){
        
        lag0 <- paste0(paste0(ind.list[[1]][[i]][[1]]$lv) , "=~", paste0(inds_label,paste0(inds),  collapse ="+"))
        lag1 <- paste0(paste0(ind.list[[1]][[i]][[1]]$lv,"lag") , "=~", paste0(inds_label,paste0(inds,"lag"), collapse ="+"))
        paste0(lag0,"\n",lag1, collapse = "")
      } else {
        paste0(ind.list[[1]][[i]][[1]]$lv ,"=~",paste0(inds, collapse="+"))
      }
    })), simplify=FALSE)
    
    } else if (scaling == "first.indicator"){
    res <- lapply(seq_along(ind.list), function(i) {
      unlist(lapply(seq_along(ind.list[[i]]), function(j) { 
        inds <- paste0(names(ind.list[[i]][[j]][[1]]$r2))
        #inds_label <- paste0(inds,"*")
        inds_label <- paste0(ind.list[[i]][[j]][[1]]$lv,"_",inds,"*")
        inds_label[1] <- ""
        if(lag){
          lag0 <- paste0(paste0(ind.list[[i]][[j]][[1]]$lv) , "=~", paste0(inds_label,paste0(inds), collapse ="+"))
          lag1 <- paste0(paste0(ind.list[[i]][[j]][[1]]$lv,"lag") , "=~", paste0(inds_label,paste0(inds,"lag"), collapse ="+"))
          paste0(lag0,"\n",lag1, collapse = "")
        } else {
          paste0(paste0(ind.list[[i]][[j]][[1]]$lv) , "=~", paste0(inds_label,paste0(inds), collapse ="+"))
        }
        
      }))
    })
    
  } else {
    
    res <- list()
    
  }
  
  
  return(res)
}