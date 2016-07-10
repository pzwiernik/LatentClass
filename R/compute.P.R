#' compute.P
#'
#' This function takes a parameter vector and computes the corresponding distribution over all variables (including the latent variable).
#' @param theta A list of length m+1, where m is the number of observed variables. The first m elements are kxr[i] matrices with the corresponding conditional probabilities. The last a element is a vector of the marginal distribution of the latent variable.
#' @details If the parameters are given in form of a vector, \code{theta.list(theta,r,n.class)} converts 
#' them to the desired format. 
#' @examples
#' theta <- list()
#' length(theta) <- 5
#' theta[[1]] <- matrix(c(0.8,1-0.7,1-0.8,0.7),2,2)
#' theta[[2]] <- matrix(c(0.8,1-0.7,1-0.8,0.7),2,2)
#' theta[[3]] <- matrix(c(0.8,1-0.7,1-0.8,0.7),2,2)
#' theta[[4]] <- matrix(c(0.8,1-0.9,1-0.8,0.9),2,2)
#' theta[[5]] <- c(1-0.7,0.7)
#' P <- compute.P(theta)
#' # to get the induced observed distribution
#' mP <- Reduce("+",P)
compute.P <- function(theta){
  # computes the fully observed distribution from the parameter theta
  # INPUT: a list of size m+1, where the first m entries are matrices of conditional distribution
  # and the last entry is the latent distribution
  # OUTPUT: a list of size k (number of hidden states), where h-th entry is the *joint* distribution
  # of X and H=h in form of a (r1,...,rm) array.
  m <- length(theta)-1
  k <- nrow(theta[[1]])
  P <- list()
  length(P) <- k
  for (h in 1:k){
    P[[h]] <- theta[[1]][h,]
    for (j in 2:m){
      P[[h]] <- outer(P[[h]],theta[[j]][h,])
    }
    P[[h]] <- theta[[length(theta)]][h]*P[[h]]
  }
  return(P)
}

