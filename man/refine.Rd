\name{refine}
\alias{refine}
\docType{package}
\title{
  Refinement step for projections
}
\description{
  \code{refine} refines a minimax design returned by \code{mMcPSO} by ensuring it has good space-filling properties on projections of the design space \eqn{X}. This refinement, in theory, does not increase the minimax distance of the design, although a slight increase may be realized in practice due to approximation error. The refined design is known as a minimax projection (miniMaxPro) design, and is further discussed in Mak and Joseph (2016). Currently only available when \eqn{X} is the unit hypercube \eqn{[0,1]^p}.
}
\usage{
  refine(mMdes,eval_pts=NA,pp_itmax=50,
          mM_tol=1e-2*(ncol(mMdes)),curscl=1,maxscl=201)
}
\arguments{
  \item{mMdes}{An \eqn{N}-by-\eqn{p} matrix for the \eqn{N}-point minimax design in \eqn{p}-dimensions.}
  \item{eval_pts}{Representative points used for approximating the continuous design space \eqn{X}.}
  \item{pp_itmax}{Maximum number of refinement iterations.}
  \item{mM_tol}{Tolerance for minimax distance increase.}
  \item{curscl,maxscl}{Redundant; still in development.}
}
\value{
  An \code{N}-by-\code{p} matrix for the minimax projection design.
}
\examples{
  \dontrun{
    #Generate a miniMaxPro design of 20 points on the unit hypercube [0,1]^8
    # (this is equivalent to runnning "des <- miniMaxPro(N=20,p=8)")
    library(foreach)
    library(randtoolbox)
    p <- 8
    mMdes <- mMcPSO(N=20,p)
    mMPdes <- refine(mMdes,eval_pts=sobol(1e6,p),mM_tol=1e-2*p)
    pairs(mMdes,xlim=c(0,1),ylim=c(0,1),pch=16) #Original minimax design from mMc-PSO
    pairs(mMPdes,xlim=c(0,1),ylim=c(0,1),pch=16) #Minimax projection design
  }
}
