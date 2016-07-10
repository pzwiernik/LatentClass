#' The log-likelihood function
#'
#' This function computes the value of the observed log-likelihood function either for a given parameter theta or for the sample distribution. 
#' @param counts The array of counts of format (r[1],...,r[m])
#' @param theta The vector of parameters in form of a list of size m+1, where the the first m entries are kxr[i] matrices of teh conditional distribution of an observed variable given the latent one. The last element is the vector of the marginal distribution of the latent variable. If not provided then the log-likelihood function at the sample distribution is computed.
#' @keywords log-likelihood function, sample distribution
#' @export 
#' @examples
#' llike(counts,theta=NULL)
llike <- function(counts,theta=NULL){
  # computes the likelihood function for the LC(r,k) model
  # INPUT:
  # counts = array of counts of format (r1,...,rm)
  # theta = a list of parameters of size m+1, where the last entry is the latent distribution
  if (is.null(theta)){
    P <- counts/sum(counts)
  }
  else if (!is.null(theta)){
    if (!is.list(theta)){
      l <- length(theta)
      r <- dim(counts)
      k <- (l+1)/(sum(r)-length(r)+1)
      theta <- theta.list(theta,r,k)
    }
    fP <- compute.P(theta)
    P <- Reduce("+",fP)  # this computes the observed distribution P(theta)
  }
  return(sum(counts*log(P)))
}

