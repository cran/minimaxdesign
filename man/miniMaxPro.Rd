\name{miniMaxPro}
\alias{miniMaxPro}
\docType{package}
\title{
  Compute minimax projection designs using clustering
}
\description{
  \code{miniMaxPro} is the main function for generating minimax projection designs on a desired design space \eqn{X}. A formal exposition of minimax projection designs and the employed algorithm can be found in Mak and Joseph (2016). Currently only available when \eqn{X} is the unit hypercube \eqn{[0,1]^p}.
}
\usage{
miniMaxPro(N,p,refine_num=1e6, refine_pts=NA, refine_itmax=50, mM_tol = 1e-02*p, ...)
}
\arguments{
  \item{N}{Number of design points desired.}
  \item{p}{Dimension of design desired.}
  \item{refine_num}{Number of representative points used to estimate \eqn{X} in the refinement step.}
  \item{refine_pts}{User-specified representative points for the refinement step. If NA, the algorithm generates these points.}
  \item{refine_itmax}{Maximum iterations for the refinement step.}
  \item{mM_tol}{Upper bound for minimax distance increase (since the refinement step may increase this distance slightly).}
  \item{...}{Parameter settings for \code{mMcPSO}.}
}
\value{
A list with two objects:
  \item{minimax}{An \code{N}-by-\code{p} matrix for the minimax design from \code{mMcPSO} before projection refinement.}
  \item{miniMaxPro}{An \code{N}-by-\code{p} matrix for the minimax projection design.}
}
\examples{
\dontrun{
#Generate a miniMaxPro design of 20 points on the unit hypercube [0,1]^8
library(foreach)
des <- miniMaxPro(N=20,p=8)
pairs(des$minimax,xlim=c(0,1),ylim=c(0,1),pch=16) #Original minimax design from mMc-PSO
pairs(des$miniMaxPro,xlim=c(0,1),ylim=c(0,1),pch=16) #Minimax projection design
}
}
