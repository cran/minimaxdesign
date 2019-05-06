\name{mMdist}
\alias{mMdist}
\docType{package}
\title{
  Computes the minimax (fill) distance of a design
}
\description{
  \code{mMdist} computes the minimax (fill) distance of a design (see Mak and Joseph (2018) for definition). Current design region options include the unit hypercube ("hypercube"), the unit simplex ("simplex"), the unit ball ("ball"), and user-defined constraints on the unit hypercube ("custom").
}
\usage{
  mMdist(D,neval=1e5,method="lattice",region="hypercube",const=NA)
}
\arguments{
  \item{D}{An \code{N}-by-\code{p} design matrix.}
  \item{neval}{Number of sample points for approximation.}
  \item{method}{Method for generating approximation points. Current options include "lattice" (lattice rule) and "sobol" (Sobol' sequence).}
  \item{region}{String for desired design region. Current options include "hypercube", "simplex", "ball", and "custom".}
  \item{const}{Function for desired constraints (inequalities) on design space.}
}
\value{
A list with two objects:
  \item{dist}{Minimax (fill) distance.}
  \item{dist.vec}{Minimax (fill) distance for each design point.}
  \item{far.pt}{A \code{p}-vector for the farthest point in design region to design.}
  \item{far.despt}{A \code{p}-vector for the design point closest to \code{far.pt}.}
}
\examples{
  \dontrun{
    #20-point minimax design on the hypercube [0,1]^2
    D <- minimax(N=20,p=2)
    plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2") #set up plot
    polygon(c(0,0,1,1),c(0,1,1,0),col="gray") #design space
    points(D,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2",pch=16) #design points
    mM <- mMdist(D)
    mM$dist #minimax (fill) distance
    lines(rbind(mM$far.pt,mM$far.despt),col="red",lty=2,lwd=2) #plot farthest point
  }
}
