#' M-step of the EM algorithm
#'
#' This internal function performs the M-step of the EM-algorithm
#' @param phat   The array of normalized counts (sample proportions) of format (r[1],...,r[m])
#' @param fphat   The list of k arrays that represent the full distribution over (r[1],...,r[m],k)
#' @keywords E-step
#' @export
#' @examples
#' M.step(phat,theta)
M.step <- function(phat,fphat){
  # the M-step of the EM algorithm
  # INPUT:
  # phat = observed sample proportions
  # fphat = full sample proportions (taken from the E-step)
  # Note: we could take only fphat because phat is easily computable but to save computing time we load both
  # OUTPUT: A new parameter vector in the standard format
  r <- dim(phat)
  m <- length(r)
  theta <- list()
  length(theta) <- m+1
  k <- length(fphat)
  for (i in 1:m){
    theta[[i]] <- matrix(0,k,r[i])
  }
  theta[[m+1]] <- rep(0,k)
  for (h in 1:k){
    theta[[m+1]][h]=sum(fphat[[h]])
    for (i in 1:m){
      theta[[i]][h,] <- apply(fphat[[h]],i,sum)/theta[[m+1]][h]
    }
  }
  return(theta)
}