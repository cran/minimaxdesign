\name{CtoB}
\alias{CtoB}
\docType{package}
\title{
  Inverse Rosenblatt transformation from the unit hypercube to the unit ball.
}
\description{
  \code{CtoB} maps points on the unit hypercube in \eqn{p}-dimensions, \eqn{C_p = [0,1]^p}, to points on the unit simplex in \eqn{p}-dimensions, \eqn{B_p}.
}
\usage{
  CtoB(D,by,num_proc)
}
\arguments{
  \item{D}{An \eqn{N}-by-\eqn{p} matrix representing \eqn{N} points in \eqn{p}-dimensions.}
  \item{by}{Step-sized used for approximating integrals.}
  \item{num_proc}{Number of processors to use.}
}
\value{
  An \eqn{N}-by-\eqn{p} matrix for the inverse-Rosenblatt mapping of \eqn{D} onto the unit ball \eqn{B_p}.
}
\examples{
  \dontrun{
    # Map the first 100 points of the Sobol' sequence in 3D
    #   onto the unit ball in 3D
    library(randtoolbox)
    des <- sobol(100,3)
    des_ball <- CtoB(des,by=1e-4,num_proc=1)

    pairs(des_ball,xlim=c(-1,1),ylim=c(-1,1),pch=16)

  }
}
