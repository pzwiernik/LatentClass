#' LC.EM function
#'
#' This function estimates parameters of the latent class model using the standard EM algorithm.
#' @param counts   The array of counts of format (r[1],...,r[m])
#' @param k   The number of latent classes fitted
#' @param tries   The number of times the EM algorithm reruns from different random starting points. The default value is tries=3
#' @param theta   The vector of parameters from which the algorithm starts. If not specified, the algorithm starts from a random point.
#' @param tol    The convergence criterion for the EM algorithm. The maximal decrease of the log-likelihood function that will terminate the algorithm.
#' @keywords EM algorithm, latent class model
#' @export
#' @examples
#' theta0 <- list()
#' length(theta0) <- 5
#' theta0[[1]] <- matrix(c(0.8,1-0.7,1-0.8,0.7),2,2)
#' theta0[[2]] <- matrix(c(0.8,1-0.7,1-0.8,0.7),2,2)
#' theta0[[3]] <- matrix(c(0.8,1-0.7,1-0.8,0.7),2,2)
#' theta0[[4]] <- matrix(c(0.8,1-0.9,1-0.8,0.9),2,2)
#' theta0[[5]] <- c(1-0.7,0.7)
#' n <- 1000
#' counts <- sample.counts(n, theta0)
#' LC.EM(counts,2)
LC.EM <- function(counts, k, tries=3, theta=NULL, tol=1e-6){
  require(gtools)
  output <- list()
  length(output) <- tries
  likes.1 <- rep(0,tries)
  r <- dim(counts)
  n <- sum(counts)
  phat <- counts/n


  st <- 1
  # if the vector of starting parameters is provided we set tries to 1.
  if (!is.null(theta)){
    st <- 0
    tries <- 1
  }

  #now we run the EM-algorithm several times starting from random starting points
  for (tr in 1:tries){
    # sample a starting point unless it was provided by the user
    if (st==1){
      theta <- sample.Theta(r,k)
    }
    l0 <- 110
    l1 <-100
    n.steps <- 0
    while (abs(l1-l0)>tol){
      l0 <- l1
      fphat <- E.step(phat,theta)
      theta <- M.step(phat,fphat)
      l1 <- llike(counts,theta)
      n.steps <- n.steps+1
    }
    likes.1[tr] <- llike(counts,theta)
    output[[tr]] <- list(parameters=theta,likelihood=likes.1[tr],iterations=n.steps)
  }
  if (tries==1){
    final <- output[[1]]
  }
  if (tries>1){
    final <- output[[which(likes.1==max(likes.1))]]
  }
  l1 <- llike(counts)
  final = c(final, LR=2*(l1-final$likelihood))
  return(final)
}


