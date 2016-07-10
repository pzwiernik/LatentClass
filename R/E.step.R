#' E-step of the EM algorithm
#'
#' This internal function performs an E-step of the EM-algorithm
#' @param phat   The array of normalized counts (sample proportions) of format (r1,...,rm)
#' @param theta   The parameter vector in form of a list of size m+1. The last element is the marginal distribution of the latent variable. The other elements of the list correspond to conditional probabilities of the observed variables given the latent one.
#' @keywords E-step
#' @export
#' @examples
#' E.step(phat,theta)
E.step <- function(phat,theta){
  # this is the E-step of the EM-algorithm
  # INPUT:
  # phat = the array of counts normalized to sum to 1 (sample proportions)
  # theta = vector of parameters
  # OUTPUT: A list of size k with the corresponding full sample proportions that 
  # maximize the conditional likelihood given the observed vector
  k <- nrow(theta[[1]])  # number of latent classes
  fP <- compute.P(theta)
  P <- Reduce("+",fP) # the induced observed distribution P(theta)
  fphat <- list()
  length(fphat) <- k
  for (h in 1:k){
    fphat[[h]] <- phat*fP[[h]]/P   # see e.g. Pachter&Sturmfels
  }
  return(fphat)
}

