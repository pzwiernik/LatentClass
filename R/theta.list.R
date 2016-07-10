#' Convert the parameter vector to a list in the internal format
#'
#' This internal function converts a parameter vector to a prefered list format 
#' @param param   The vector of parameters of the model. It represents elements of each transition matrix apart from the last column.
#' @param r The dimension vector of the data (r[1],...,r[m])
#' @param n.class The number of hidden classes
#' @keywords theta
#' @export
#' @examples
#' theta.list(param,r,class=2)
theta.list <- function(param,r,n.class=2){
  # the length of param should be k*(r1-1)+...+k*(rm-1)+k-1=k*(R-m+1)-1
  theta <- list()
  index <- 1
  m <- length(r)
  for (i in 1:m){
    index1 <- index+n.class*(r[i]-1)-1
    theta[[i]] <- matrix(param[index:index1],n.class,r[i]-1)
    row.sum <- apply(theta[[i]],1,sum)
    theta[[i]] <- cbind(theta[[i]],1-row.sum)
    index <- index1+1
  }
  index1 <- index+n.class-2
  theta[[m+1]] <- c(param[index:index1],1-sum(param[index:index1]))
  return(theta)
}

  