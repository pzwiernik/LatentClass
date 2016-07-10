219.15+175+8.50+8.30+8+8+18.30+8+19.50+1.30+12+6.50+33+11.30#' Sample a parameter vector
#'
#' This function samples a parameter vector for a latent class model. Each row from each conditional distribution matrix is sampled from the uniform Dirichlet distribution. 
#' @param r A vector  (r[1],...,r[m]) of integers, where m is the number of observed variables and r[i] is teh statespace of the i-th variable.
#' @param k The number of latent classes fitted
#' @keywords sample parameters
#' @export
#' @examples
#' r <- c(2,2,2,2)
#' k <- 3
#' sample.Theta(r,k)
sample.Theta <- function(r,k){
  # This function gives a random vector of parameters for a LC(r,k) model
  # Input: k - number of latent classes, r=(r1,...,rm) number of observed states.
  # Output:  a list of size (m+1)
  theta <- list()
  m <- length(r)
  length(theta) <- m+1
  for (i in 1:m){
    theta[[i]] <- rdirichlet(k, rep(1,r[i]))
  }
  theta[[m+1]] <- (k:1)*2/(k*(k+1))
  return(theta)
}