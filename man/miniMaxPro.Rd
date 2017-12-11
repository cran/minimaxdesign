\name{miniMaxPro}
\alias{miniMaxPro}
\docType{package}
\title{
  Compute minimax projection designs using clustering
}
\description{
  \code{miniMaxPro} is the main function for generating minimax projection designs on a desired design space \eqn{X}. A formal exposition of minimax projection designs and the employed algorithm can be found in Mak and Joseph (2017). Currently only available when \eqn{X} is the unit hypercube \eqn{[0,1]^p}.
}
\usage{
miniMaxPro(N,p,mMdes=NA, mM_tol=1e-3*p,
            refine_num=1e5, refine_pts=NA, refine_itmax=100, ...)
}
\arguments{
  \item{N}{Number of design points desired.}
  \item{p}{Dimension of design desired.}
  \item{mMdes}{Minimax design from \code{mMcPSO()}.}
  \item{refine_num}{Number of points used to estimate \eqn{X} in the refinement step.}
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
des <- miniMaxPro(N=40,p=8,refine_itmax=100)
pairs(des$minimax,xlim=c(0,1),ylim=c(0,1),pch=16)
pairs(des$miniMaxPro,xlim=c(0,1),ylim=c(0,1),pch=16)

#Use the minimax design from mMc-PSO to warm start miniMaxPro at a new setting
mMdes <- mMcPSO(N=40,p=8)
newdes <- miniMaxPro(N=20,p=8,mMdes=mMdes,refine_itmax=100)
pairs(newdes$miniMaxPro,xlim=c(0,1),ylim=c(0,1),pch=16)

}
}
