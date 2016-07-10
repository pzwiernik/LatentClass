#' M-step of the EM algorithm
#'
#' This internal function converts a parameter vector to a prefered list format 
#' @param phat   The array of normalized counts (sample proportions) of format (r[1],...,r[m])
#' @param fphat   The list of k arrays that represent the full distribution over (r[1],...,r[m],k)
#' @keywords E-step
#' @export
#' @examples
#' theta.list(param,r,class=2)
theta.unlist <- function(theta){
  # the length of param should be k*(r1-1)+...+k*(rm-1)+k-1=k*(R-m+1)-1
  m <- length(theta)-1
  param <- c()
  for (i in 1:m){
    ri <- ncol(theta[[i]])
    param <- c(param,c(theta[[i]][,1:(ri-1)])) 
  }
  k <- length(theta[[m+1]])
  param <- c(param,c(theta[[m+1]][1:(k-1)])) 
  return(param)
}

