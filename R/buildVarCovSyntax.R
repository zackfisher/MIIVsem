#' build the variance covariance point estimates model syntax
#'@keywords internal
buildVarCovSyntax <- function(pt){
  
  # model syntax
  var.cov.model <- unlist(apply(pt, 1, function(r){
    
    if (r["op"] %in% c("~","=~")) {
      
      paste0(r["lhs"],r["op"],r["ustart"],"*",r["rhs"])
      
    } else if (r["op"] == "~1"){
      
      paste0(r["lhs"],"~",r["ustart"],"*", "1")
      
    } else if (r["op"] == "~~" & r["label"] != ""){
      
      paste0(r["lhs"],"~~",r["label"],"*", r["rhs"])
      
    } else if (r["op"] == "~~" & r["label"] == ""){
      
      paste0(r["lhs"],"~~", r["rhs"])
      
    }
    
  }))
  
  return(var.cov.model)
  
}
