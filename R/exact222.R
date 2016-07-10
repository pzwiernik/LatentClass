#' Exact maximum likelihood for three binary variables
#'
#' This function computes the exact maximum likelihood for the model with three observed binary variables and two latent classes. 
#' @param counts   A (2,2,2) array with counts. 
#'        tol   tolerance with respect to which model membership is checked (default tol=1e-8)
#' @keywords exact MLE
#' @export
#' @examples
#' n <-  c(12, 23, 34, 45, 56, 67, 78, 90)
#' res <- exact222(n)
#' res$min.deviance
#' res$fitted.values

exact222 <- function(counts,tol=1e-8){

  states <-  factor(0:1)
  data <- expand.grid(x1 = states, x2 = states, x3 = states)
  data <- data.frame(Freq = counts, data)
  options(contrasts = c('contr.treatment', 'contr.poly'))
  mdevs <- rep(Inf,30)
  fitted <- list()
  length(fitted) <- 30
  
  counter <- 0
  Z <- model.matrix(~ x1*x2*x3, data = data)[,]
  mo <- glm(Freq ~ 0 + Z, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z1 <- cbind(Z[,1:6], Z[,8])
  mo <- glm(Freq ~ 0 + Z1, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z2 <- cbind(Z[,1:6], Z[,7] - Z[,8])
  mo <- glm(Freq ~ 0 + Z2, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z3 <- cbind(Z[,1:5], Z[,7],Z[,8])
  mo <- glm(Freq ~ 0 + Z3, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z4 <- cbind(Z[,1:5], Z[,6] - Z[,8], Z[,7])
  mo <- glm(Freq ~ 0 + Z4, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z5 <- cbind(Z[,1:4], Z[,6:8])
  mo <- glm(Freq ~ 0 + Z5, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z6 <- cbind(Z[,1:4], Z[,5] - Z[,8], Z[,6:7])
  mo <- glm(Freq ~ 0 + Z6, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  
  Z12 <- cbind(Z[,1:6])
  mo <- glm(Freq ~ 0 + Z12, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z13 <- cbind(Z[,1:5],Z[,8])
  mo <- glm(Freq ~ 0 + Z13, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z14 <- cbind(Z[,1:5],Z[,6]-Z[,8])
  mo <- glm(Freq ~ 0 + Z14, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z15 <- cbind(Z[,1:4],Z[,6],Z[,8])
  mo <- glm(Freq ~ 0 + Z15, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z16 <- cbind(Z[,1:4],Z[,5]-Z[,8],Z[,6])
  mo <- glm(Freq ~ 0 + Z16, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z23 <- cbind(Z[,1:5],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z23, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z24 <- cbind(Z[,1:5],Z[,6]-Z[,8],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z24, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z25 <- cbind(Z[,1:4],Z[,6],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z25, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z26 <- cbind(Z[,1:4],Z[,5]-Z[,8],Z[,6],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z26, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z34 <- cbind(Z[,1:5],Z[,7])
  mo <- glm(Freq ~ 0 + Z34, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z35 <- cbind(Z[,1:4],Z[,7:8])
  mo <- glm(Freq ~ 0 + Z35, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z36 <- cbind(Z[,1:4],Z[,5]-Z[,8],Z[,7])
  mo <- glm(Freq ~ 0 + Z36, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z45 <- cbind(Z[,1:4],Z[,6]-Z[,8],Z[,7])
  mo <- glm(Freq ~ 0 + Z45, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z46 <- cbind(Z[,1:4],Z[,5]-Z[,8],Z[,6]-Z[,8],Z[,7])
  mo <- glm(Freq ~ 0 + Z46, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z56 <- cbind(Z[,1:4],Z[,6:7])
  mo <- glm(Freq ~ 0 + Z56, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z135 <- cbind(Z[,1:4],Z[,8])
  mo <- glm(Freq ~ 0 + Z135, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z136 <- cbind(Z[,1:4],Z[,5]-Z[,8])
  mo <- glm(Freq ~ 0 + Z136, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z145 <- cbind(Z[,1:4],Z[,6]-Z[,8])
  mo <- glm(Freq ~ 0 + Z145, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z235 <- cbind(Z[,1:4],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z235, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z146 <- cbind(Z[,1:4],Z[,5]-Z[,8],Z[,6]-Z[,8])
  mo <- glm(Freq ~ 0 + Z146, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z236 <- cbind(Z[,1:4],Z[,5]-Z[,8],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z236, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z245 <- cbind(Z[,1:4],Z[,6]-Z[,8],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z245, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }
  
  Z246 <- cbind(Z[,1:4],Z[,5]-Z[,8],Z[,6]-Z[,8],Z[,7]-Z[,8])
  mo <- glm(Freq ~ 0 + Z246, poisson, data= data)
  counter <- counter+1
  if (test222(array(mo$fitted.values,c(2,2,2)),tol)==1) {
    mdevs[counter] <- mo$deviance
    fitted[[counter]] <- mo$fitted.values
  }

  strata <- list(c(),1,2,3,4,5,6,c(1,2),c(1,3),c(1,4),c(1,5),c(1,6),c(2,3),c(2,4),
                 c(2,5),c(2,6),c(3,4),c(3,5),c(3,6),c(4,5),c(4,6),c(5,6),c(1,3,5),
                 c(1,3,6),c(1,4,5),c(2,3,5),c(1,4,6),c(2,3,6),c(2,4,5),c(2,4,6))
  i.max <- which(mdevs==min(mdevs))
  return(list(min.deviance=mdevs[i.max],fitted.values=fitted[[i.max]],
              deviances=mdevs,m.strata=strata[[i.max]],strata=strata))
}