\name{CtoB2}
\alias{CtoB2}
\docType{package}
\title{
  Inverse Rosenblatt transformation from the unit hypercube to the unit ball in two-dimensions.
}
\description{
  \code{CtoB2} maps points on the unit hypercube in \eqn{2}-dimensions, \eqn{C_2 = [0,1]^2}, to points on the unit ball in \eqn{2}-dimensions, \eqn{B_2}.
}
\usage{
  CtoB2(D,by,num_proc)
}
\arguments{
  \item{D}{An \eqn{N}-by-\eqn{p} matrix representing \eqn{N} points in \eqn{2}-dimensions.}
  \item{by}{Step-sized used for approximating integrals.}
  \item{num_proc}{Number of processors to use.}
}
\value{
  An \eqn{N}-by-\eqn{p} matrix for the inverse-Rosenblatt mapping of \eqn{D} onto the unit ball \eqn{B_2}.
}
\examples{
  \dontrun{
    # Map the first 100 points of the Sobol' sequence in 2D
    #   onto the unit ball in 2D
    library(randtoolbox)
    des <- sobol(100,2)
    des_ball <- CtoB2(des,by=1e-4,num_proc=1)

    plot(des_ball,xlim=c(-1,1),ylim=c(-1,1),pch=16)

  }
}
