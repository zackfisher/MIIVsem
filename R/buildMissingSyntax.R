#' build the missing data model syntax
#'@keywords internal
buildMissingSyntax <- function(pt){
  
  missing.model <- unlist(apply(pt, 1, function(r){
    
    # operate on the results of estimatesTable, assuming:
    
    # 1. if mlabel is not NA:
    
    if (!is.na(r["mlabel"])){
      
      # 1.a. mlabel is numeric and there is a fixed coefficient, 
      #      in which case the parameter should be fixed
      
      if (is.numeric(utils::type.convert(r["mlabel"], as.is=TRUE))){
        
        if (r["op"] %in% "~1"){
          paste(r["lhs"], "~",r["mlabel"],"*1\n",collapse="")
        } else {
          paste(r["lhs"],r["op"],r["mlabel"],"*",r["rhs"],"\n",collapse="")
        }
        
        # 1.b. there is a label, in which case the label should be
        #      assigned to the parameter. we should also assign
        #      starting values here if applicable. 
        
      } else {
        
        if (r["op"] %in% "~1" & as.numeric(r["free"]) == 0){
          
          paste(r["lhs"], "~",r["ustart"],"*1\n",collapse="")
          
        } else if (r["op"] %in% "~1" & as.numeric(r["free"]) != 0){
          
          paste(r["lhs"], "~",r["mlabel"],"*1\n",collapse="")
          paste(r["lhs"], "~","start(",r["ustart"],")*1\n",collapse="")
          
        } else if (r["op"] %in% c("~", "~~","=~") & as.numeric(r["free"]) == 0){
          
          paste(r["lhs"],r["op"],r["ustart"],"*",r["rhs"],"\n",collapse="")
          
        } else if (r["op"] %in% c("~", "~~","=~") & as.numeric(r["free"]) != 0){
          
          paste(r["lhs"],r["op"],r["mlabel"],"*",r["rhs"],"\n",collapse="")
          paste(r["lhs"],r["op"],"start(",r["ustart"],")*",r["rhs"],"\n",collapse="")
        }
      } 
      
      # 2. mlabel is NA:
    } else {
      if(r["op"] %in% c("~~","~","=~","~1")){
        if (as.numeric(r["free"]) == 0 ){
          if (r["op"] %in% "~1"){
            paste(r["lhs"], "~",r["ustart"],"*1\n",collapse="")
          } else {
            paste(r["lhs"],r["op"],r["ustart"],"*",r["rhs"],"\n",collapse="")
          }
        } else if (as.numeric(r["free"] != 0 )){
          if (r["op"] %in% "~1"){
            paste(r["lhs"], "~start(",r["ustart"],")*1\n",collapse="")
          } else {
            paste(r["lhs"],r["op"],"start(",r["ustart"],")*",r["rhs"],"\n",collapse="")
          }
        }
      }
    }
  }))
  
  return(missing.model)
  
}
