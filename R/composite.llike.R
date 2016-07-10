#' Compute the composite likelihood
#'
#' This internal function computes the composite likelihood for a given latent class model 
#' and a collection of margins.
#' @param dat The data array of the format \code{(r[1],...,r[m])}
#' @param param The vector of parameters of the model. It represents elements of each transition matrix apart from the last column.
#' @param n.class The number of hidden classes
#' @param margin.size The size of the margins that are used to compute the composite lieklihood.
composite.llike <- function(dat,param,n.class=2,margin.size=3){
  # computes minus the composite likelihood
  r <- dim(dat)
  theta <- theta.list(param,r,n.class)
  fP <- compute.P(theta)
  P <- Reduce("+",fP)  # this computes the observed distribution P(theta)
  margins <- combn(length(r),margin.size)
  n.margs <- ncol(margins)
  mP <- list()
  mU <- list()
  mL <- list()
  length(mP) <- length(mU) <- length(mL) <- n.margs
  for (i in 1:n.margs){
    mP[[i]] <- apply(P,margins[,i],sum)
    mU[[i]] <- apply(dat,margins[,i],sum)
    mL[[i]] <- sum(mU[[i]]*log(mP[[i]]))
  }
  return(Reduce("+",mL))
}

