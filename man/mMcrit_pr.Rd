\name{mMcrit_pr}
\alias{mMcrit_pr}
\docType{package}
\title{
  Compute minimax distances in projections
}
\description{
  \code{mMcrit_pr} computes the maximum minimax distance of a design \eqn{D} when the desired design space \eqn{X}is projected onto all possible subspaces of dimension \eqn{k} (see Mak and Joseph, 2016 for a formal definition). Currently only available when \eqn{X} is the unit hypercube \eqn{[0,1]^p}.
}
\usage{
mMcrit_pr(D,k,eval_num)
}
\arguments{
  \item{D}{An \code{N}-by-\code{p} matrix representing an \\code{N}-point design in \code{p}-dimensions.}
  \item{k}{Projection dimension.}
  \item{eval_num}{Number of representative points used to estimate projected subspaces.}
}
\value{
The maximum estimated minimax distance of the design \code{D} over all possible projections of the design space \code{X} with dimension \eqn{k}.
}
\examples{
\dontrun{
#Estimate the projected minimax distance in 2 dimensions
# of a 20-point Sobol' sequence design in 8 dimensions.
library(randtoolbox)
D <- sobol(20,8)
mMcrit_pr(D,2,1e7) #1e7 points from Sobol' sequence used in estimate
}
}
