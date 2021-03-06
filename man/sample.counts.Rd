% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/sample.counts.R
\name{sample.counts}
\alias{sample.counts}
\title{Sample data}
\usage{
sample.counts(n, theta)
}
\arguments{
\item{n}{The sample size}

\item{theta}{The parameter vector in form of a list...}
}
\description{
This samples data from a distribution in the latent class model
}
\examples{
theta0 <- list()
length(theta0) <- 5
theta0[[1]] <- matrix(c(0.8,0.8,1-0.8,1-0.8),2,2) # rank 1
theta0[[2]] <- matrix(c(0.8,0.8,1-0.8,1-0.8),2,2) # rank 1
theta0[[3]] <- matrix(c(0.8,0.8,1-0.8,1-0.8),2,2) # rank 1
theta0[[4]] <- matrix(c(0.8,0.4,1-0.8,1-0.4),2,2)
theta0[[5]] <- c(1-0.3,0.3)
n <- 10000
counts <- sample.counts(n, theta0)
}
\keyword{cats}

