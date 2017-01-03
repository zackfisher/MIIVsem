#' @keywords internal
calcStage1Fitted <- function(Z_i, V_i) {
  
  Zhat_i <- list()
  
  for(i in 1:length(Z_i)) {
    stage1 <- lm.fit(as.matrix(V_i[[i]]), as.matrix(Z_i[[i]]))
    Zhat_i[[i]] <- as.matrix(stage1$fitted.values)
    #Zhat_i[[i]] <- V_i[[i]] %*% solve(crossprod(V_i[[i]]), crossprod(V_i[[i]], Z_i[[i]]))
  }
  return(Zhat_i)
}