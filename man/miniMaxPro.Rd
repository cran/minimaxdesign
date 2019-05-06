\name{miniMaxPro}
\alias{miniMaxPro}
\docType{package}
\title{
  Compute minimax projection designs using clustering on constrained design regions
}
\description{
  \code{miniMaxPro}.
}
\usage{
miniMaxPro(N,p,mMdes=NA, mMtol=1e-3*p,
          neval=1e5, itmax_refine=100, ...)
}
\arguments{
  \item{N}{Number of design points.}
  \item{p}{Dimension of design region.}
  \item{mMdes}{Minimax design from \code{minimax()}.}
  \item{mMtol}{Tolerance for fill distance increase (actual increase may be slightly higher).}
  \item{neval}{Number of sample points for refinement.}
  \item{itmax_refine}{Maximum number of iterations for refinement.}
  \item{...}{Parameters for \code{minimax()}.}
}
\value{
A list with two objects:
  \item{minimax}{An \code{N}-by-\code{p} matrix for the minimax design.}
  \item{miniMaxPro}{An \code{N}-by-\code{p} matrix for the minimax projection design.}
}
\examples{
\dontrun{
#30-point miniMaxPro design on the hypercube [0,1]^6
D <- minimax(N=30,p=6)
D <- miniMaxPro(N=30,p=6,mMdes=D)
mMdist(D$minimax)$dist
mMdist(D$miniMaxPro)$dist #slightly higher fill distance
pairs(D$minimax,xlim=c(0,1),ylim=c(0,1),pch=16)
pairs(D$miniMaxPro,xlim=c(0,1),ylim=c(0,1),pch=16) #... but better projections

}
}
