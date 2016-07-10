#' constraints.theta function
#'
#' This auxiliary function outputs for any latent class model the corresponding matrix that encodes paramameter constraints. 
#' The output is a list \code{A,b}, where \code{A} is a matrix and \code{b} is a vector such that \code{A\%*\% theta >= b} is the set of the inequalities defining parameter constraints.
#' @param r   The vector of dimensions of the format \code{(r[1],...,r[m])}
#' @param n.class   The number of latent classes
#' @details In addition, for identifiability issues, we require that teh latent probabilities are sorted
#' from the largest to the smallest.
#' @examples
#' constraints.theta(c(2,2,2,2),2)

constraints.theta <- function(r,n.class){
  m <- length(r)
  lpar <- n.class*(sum(r)-m+1)-1 #this is the length of the parameter vector
  A <- diag(lpar)
  b <- rep(0,lpar)
  
  index <- 1
  for (i in 1:m){
    for (j in 1:n.class){
      vect <- rep(0,lpar) 
      index1 <- index+r[i]-2
      vect[index:index1] <- rep(-1,r[i]-1)
      index <- index1+1
      A <- rbind(A,vect)
      b <- c(b,-1)
    }
  }
  vect <- rep(0,lpar) 
  vect[index:(index+(n.class-2))] <- rep(-1,n.class-1)
  A <- rbind(A,vect)
  b <- c(b,-1) 
  if (n.class==2){
    vect <- rep(0,lpar)
    vect[lpar] <- 2
    A <-rbind(A,vect)
    b <- c(b,1)
  }
  if (n.class>2){
    for (i in 1:(n.class-2)){
      vect <- rep(0,lpar) 
      vect[(n.class*(sum(r)-m)+i):(n.class*(sum(r)-m)+i+1)] <- c(1,-1) 
      A <-rbind(A,vect)
      b <- c(b,0)
    }
    vect <- rep(0,lpar)
    vect[(n.class*(sum(r)-m)+1):(lpar-1)] <- rep(1,k-2) 
    vect[lpar] <- 2
    A <-rbind(A,vect)
    b <- c(b,1)
  }
 return(list(A,b))
}