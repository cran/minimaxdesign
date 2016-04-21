\name{CtoA}
\alias{CtoA}
\docType{package}
\title{
  Inverse Rosenblatt transformation from the unit hypercube to the unit simplex
}
\description{
  \code{CtoA} maps points on the unit hypercube in \eqn{p}-dimensions, \eqn{C_p = [0,1]^p}, to points on the unit simplex in \eqn{p}-dimensions, \eqn{A_p}.
}
\usage{
  CtoA(D,by,num_proc)
}
\arguments{
  \item{D}{An \eqn{N}-by-\eqn{p} matrix representing \eqn{N} points in \eqn{p}-dimensions.}
  \item{by}{Step-sized used for approximating integrals.}
  \item{num_proc}{Number of processors to use.}
}
\value{
  An \eqn{N}-by-\eqn{p} matrix for the inverse-Rosenblatt mapping of \eqn{D} onto the unit simplex \eqn{A_p}.
}
\examples{
  \dontrun{
    # Map the first 100 points of the Sobol' sequence in 2D
    #   onto the unit simplex in 2D
    library(randtoolbox)
    des <- sobol(100,2)
    des_simp <- CtoA(des,by=1e-4,num_proc=1)

    plot(des_simp,xlim=c(0,1),ylim=c(0,1),pch=16)

  }
}
