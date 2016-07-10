#' LC.CL function
#'
#' This function estimates parameters of the latent class model using the composite likelihood.
#' As a subroutine it uses \code{constrOptim()}.
#' @param counts   The array of counts of format \code{(r[1],...,r[m])}
#' @param k   The number of latent classes fitted
#' @param tries   The number of times the algorithm reruns from different random starting points. The default value is tries=3
#' @param theta   The vector of parameters from which the algorithm starts. If not specified, the algorithm starts from a random point. 
#' @param tol    The convergence criterion for the EM algorithm. The maximal decrease of the log-likelihood function that will terminate the algorithm.
#' @param margin.size The size of the margins that will be included in the construction of the composite likelihood. 
#' @details In this first version there is no way to specify the margins that will be taken 
#' into account in the construction of the composite likelihood function. Instead, for a given 
#' \code{margin.size} all margins of this size are taken into account. The defult is \code{3} because the
#' latent class model is identified from all the triples.
#' @export 
#' @examples
#' theta0 <- list()
#' theta0[[1]] <- matrix(c(0.9,0.2,0.1,0.8),2,2)
#' theta0[[2]] <- matrix(c(0.7,0.2,0.3,0.8),2,2)
#' theta0[[3]] <- matrix(c(0.9,0.1,0.1,0.9),2,2)
#' theta0[[4]] <- matrix(c(0.7,0.1,0.3,0.9),2,2)
#' theta0[[5]] <- matrix(c(0.7,0.1,0.3,0.9),2,2)
#' theta0[[6]] <- c(0.3,0.7)
#' data <- sample.counts(500,theta0)
#' LC.CL(data,2,tries=3,tol=1e-8,margin.size=3)
LC.CL <- function(counts, k, tries=3, theta=NULL, tol=1e-8, margin.size=3){
  output <- list()
  length(output) <- tries  
  likes.1 <- rep(0,tries)
  r <- dim(counts)

  st <- 1
  # if the vector of starting parameters is provided we set tries to 1. 
  if (!is.null(theta)){
    st <- 0
    tries <- 1
  }
  constr <- constraints.theta(r,k)
  #now we run the optimization procedure several times starting from random starting points
  for (tr in 1:tries){
    # sample a starting point unless it was provided by the user
    if (st==1){
      theta <- sample.Theta(r,k)    
    }
    param <- theta.unlist(theta)
    cres <- constrOptim(param,composite.llike,grad=NULL,ui=constr[[1]],ci=constr[[2]],dat=data,margin.size=margin.size,control=list(reltol=tol,fnscale=-1))
#    cres <- constrOptim(cres$par,composite.llike,grad=NULL,ui=constr[[1]],ci=constr[[2]],dat=data,margin.size=margin.size)
    likes.1[tr] <- cres$value  
    output[[tr]] <- list(parameters=cres$par,likelihood=likes.1[tr])
  }
  if (tries==1){
    final <- output[[1]]
  }
  if (tries>1){
    final <- output[[which(likes.1==max(likes.1))]]
  }
  return(final)
}