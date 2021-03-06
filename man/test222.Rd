% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/test222.R
\name{test222}
\alias{test222}
\title{Internal function to check if a 2x2x2 distribution lies in LC(2) model}
\usage{
test222(p, tol = 1e-08)
}
\arguments{
\item{counts}{A (2,2,2) array with counts.}
}
\description{
This function check if a 2x2x2 distribution lies in LC(2) model.
}
\examples{
r <- c(2,2,2,2)
k <- 3
sample.Theta(r,k)
}
\keyword{MLE}
\keyword{exact}

