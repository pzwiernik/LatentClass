% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/LC.MLE.R
\name{LC.MLE}
\alias{LC.MLE}
\title{LC.MLE function}
\usage{
LC.MLE(counts, k, tries = 3, theta = NULL, tol = 1e-08, mode = "optim")
}
\arguments{
\item{counts}{The array of counts of format \code{(r[1],...,r[m])}}

\item{k}{The number of latent classes fitted}

\item{tries}{The number of times the algorithm reruns from different random starting points. The default value is tries=3}

\item{theta}{The vector of parameters from which the algorithm starts. If not specified, the algorithm starts from a random point.}

\item{tol}{The convergence criterion for the EM algorithm. The maximal decrease of the log-likelihood function that will terminate the algorithm.}

\item{mode}{The algorithm used to compute the MLE. Write \code{'EM'} to use the EM algorithm and \code{'optim'}
to use the standard R constrained optimization routine \code{constrOptim()}.}
}
\description{
This function estimates parameters of the latent class model using the maximum likelihood method.
As a subroutine it uses \code{constrOptim()}.
}
\details{
This function works as an alternative to \code{LC.EM()} that uses the EM algorithm.
}
\examples{
theta0 <- list()
theta0[[1]] <- matrix(c(0.9,0.2,0.1,0.8),2,2)
theta0[[2]] <- matrix(c(0.7,0.2,0.3,0.8),2,2)
theta0[[3]] <- matrix(c(0.9,0.1,0.1,0.9),2,2)
theta0[[4]] <- matrix(c(0.7,0.1,0.3,0.9),2,2)
theta0[[5]] <- matrix(c(0.7,0.1,0.3,0.9),2,2)
theta0[[6]] <- c(0.3,0.7)
data <- sample.counts(500,theta0)
LC.MLE(data,2,tries=3,tol=1e-8,margin.size=3)
}

