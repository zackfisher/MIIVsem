#' @keywords internal
model.matrix.miive <- function(object, which = "x", ... ){
   
  result <- matrix(NA, 0, 0)
  mmRowNames <- NULL
  mmColNames <- NULL
  
  for(i in 1:length(object$eq)) {
    
    mmi <- model.matrix(object$eq[[i]], which = which)
    
    
    result <- rbind(
         cbind( result, matrix( 0, nrow( result ), ncol( mmi ) ) ),
         cbind( matrix( 0, nrow( mmi ), ncol( result ) ), mmi ) )
    
      mmRowNames <- c( mmRowNames,
         paste( object$eq[[ i ]]$eqnLabel, "_", rownames( mmi ), sep = "" ) )
      for( j in 1:ncol( mmi ) ){
         cName <- colnames( mmi )[ j ]
         if( object$panelLike && cName != "(Intercept)" ){
            mmColNames <- c( mmColNames, cName )
         } else {
            mmColNames <- c( mmColNames,
               paste( object$eq[[ i ]]$eqnLabel, "_", cName, sep = "" ) )
         }
      }
   }
   rownames( result ) <- mmRowNames
   colnames( result ) <- mmColNames
   return( result )
}

## return model matrix of a single equation
model.matrix.miive.equation <- function(object, which = "x", ...){
   
  # Throw a warning if wrong which is used.
  if(!which %in% c( "Z", "Zhat", "V")) {
      stop("MIIVsem: Argument 'which' must be  \"Z\", \"Zhat\", or \"V\"" )
  } 
  
  if(which == "ZHat"){
    
    Z   <- model.matrix(object, which = "Z")
    V   <- model.matrix(object, which = "V")
    res <- residuals( object )
      
    if(sum(!is.na(res)) != nrow(Z)) {
      stop("MIIVsem internal error: # nrow(residuals) != nrow(Z)")
    }else if (nrow(Z) != nrow(V)) {
      stop("MIIVsem internal error: # nrow(V) != nrow(Z)")
    }
    
    result <- V %*% solve(crossprod(V), crossprod(V, X))
  
  } else {
      
    if (!is.null(object[[which]])) {
      
      result <- object[[which]]
      
    } else if(!is.null(model.frame(object))) {
    
      if(which == "Z") {
        
        result <- model.matrix(object$terms_zy_i, data = model.frame(object))
        
      } else {
            
        result <- model.matrix(object$terms_v_i, data = object$modelInst )
         
      }
      
      attrAssign <- attributes(result)$assign
      result <- result[!is.na(residuals(object)), , drop = FALSE]
      attributes(result)$assign <- attrAssign
         
      } else {
         # Nothing here.
      }
   }

   return( result )
}
