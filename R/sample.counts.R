#' Sample data 
#'
#' This samples data from a distribution in the latent class model
#' @param n The sample size
#' @param theta The parameter vector in form of a list...
#' @keywords cats
#' @export
#' @examples
#' theta0 <- list()
#' length(theta0) <- 5
#' theta0[[1]] <- matrix(c(0.8,0.8,1-0.8,1-0.8),2,2) # rank 1
#' theta0[[2]] <- matrix(c(0.8,0.8,1-0.8,1-0.8),2,2) # rank 1
#' theta0[[3]] <- matrix(c(0.8,0.8,1-0.8,1-0.8),2,2) # rank 1
#' theta0[[4]] <- matrix(c(0.8,0.4,1-0.8,1-0.4),2,2) 
#' theta0[[5]] <- c(1-0.3,0.3)
#' n <- 10000
#' counts <- sample.counts(n, theta0)
sample.counts <- function(n,theta){
  # sample data from the given distribution P(theta)
  fP <- compute.P(theta)
  P <- Reduce("+",fP)  # this computes the observed distribution P(theta)
  m <- length(theta)-1
  r <- rep(0,m)
  for (i in 1:m){
    r[i] <- ncol(theta[[i]])
  }
  nx <- prod(r)
  dat <- sample(1:nx,n,replace=TRUE,c(P)) 
  counts <- tabulate(dat)
  dim(counts) <- r
  return(counts)
}

