\name{mMcrit}
\alias{mMcrit}
\docType{package}
\title{
  Compute minimax distances
}
\description{
  \code{mMcrit} computes the minimax distance of a design \eqn{D} to the desired design space \eqn{X} (see Mak and Joseph, 2016 for a formal definition). A large set of representative points is used to estimate \eqn{X} when it is a continuous space.
}
\usage{
mMcrit(D,eval_pts)
}
\value{
A list with three objects:
  \item{mM_dist}{The estimated minimax distance of design \eqn{D}.}
  \item{far_pt}{The farthest point in the design space \eqn{X} to the design \eqn{D}.}
  \item{des_pt}{The closest design point to \code{far_pt}.}
}
\arguments{
  \item{D}{An \code{N}-by-\code{p} matrix representing an \\code{N}-point design in \code{p}-dimensions.}
  \item{eval_pts}{A large matrix of representative points for estimating the design space \eqn{X}. This can be obtained by uniform random sampling on \eqn{X}, or using low-discrepancy sequences.}
}
\examples{

\dontrun{
#Estimate the minimax distance of a 20-point
# Sobol' sequence design in 8 dimensions.
library(randtoolbox)
D <- sobol(20,8)
mMcrit(D,sobol(1e7,8)) #1e7 points from Sobol' sequence used in estimate
}
}
