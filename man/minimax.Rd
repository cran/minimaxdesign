\name{minimax}
\alias{minimax}
\docType{package}
\title{
  Compute minimax designs using clustering on constrained design regions
}
\description{
  \code{minimax} computes minimax designs using the minimax clustering algorithm in Mak and Joseph (2018). Current design region options include the unit hypercube ("hypercube"), the unit simplex ("simplex"), the unit ball ("ball"), and user-defined constraints on the unit hypercube ("custom").
}
\usage{
minimax(N,p,q=10,region="hypercube",ini=NA,const=NA,
        params_pso=list(w=0.72,c1=1.49,c2=1.49),
        npart=5,nclust=1e5,neval=nclust,
        itmax_pso=50,itmax_pp=100,itmax_inn=1e4,jit=0.1/sqrt(N))
}
\arguments{
  \item{N}{Number of design points.}
  \item{p}{Dimension of design region.}
  \item{q}{Power parameter for approximating the minimax criterion (see paper for details). Larger values of q give a better approximation, but may cause numerical instability.}
  \item{region}{String for desired design region. Current options include "hypercube", "simplex", "ball", and "custom".}
  \item{ini}{String for initial design. Current choices include NA (automatic choice), "sobol" (Sobol' sequence), and "ff" (fractional factorial).}
  \item{const}{Function for desired constraints (inequalities) on design space. See examples for implementation.}
  \item{params_pso}{Particle swarm optimization parameters (particle momentum (w), local-best velocity (c1) and global-best velocity (c2)).}
  \item{npart}{Number of particles for particle swarm optimization.}
  \item{nclust,neval}{Number of sample points for minimax clustering and post-processing.}
  \item{itmax_pso,itmax_pp,itmax_inn}{Maximum number of iterations for minimax clustering, post-processing and inner optimization.}
  \item{jit}{Jitter radius for post-processing.}
}
\value{
  An \code{N}-by-\code{p} matrix for the minimax design.
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

#20-point minimax design on the simplex A_2
D <- minimax(N=20,p=2,region="simplex")
plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2") #set up plot
polygon(c(0,0,1),c(0,1,1),col="gray") #design space
points(D,pch=16) #design points
mM <- mMdist(D,region="simplex")
mM$dist #minimax (fill) distance
lines(rbind(mM$far.pt,mM$far.despt),col="red",lty=2,lwd=2) #plot farthest point

#20-point minimax design on the ball B_2
library(plotrix)
D <- minimax(N=20,p=2,region="ball")
plot(NULL,xlim=c(-1,1),ylim=c(-1,1),xlab="x1",ylab="x2") #set up plot
draw.circle(0,0,1,col="gray") #design space
points(D,pch=16) #design points
mM <- mMdist(D,region="ball")
mM$dist #minimax (fill) distance
lines(rbind(mM$far.pt,mM$far.despt),col="red",lty=2,lwd=2) #plot farthest point

#20-point minimax design on custom design space (inequalities on [0,1]^2)
ineqs <- function(xx){ #user-defined inequalities
  bool.vec <- rep(TRUE,length(xx))
  bool.vec[1] <- (xx[2]<=2-2*xx[1]) #inequality 1: x2 <= 2 - 2*x1
  bool.vec[2] <- (xx[1]>=xx[2]) #inequality 2: x1 >= x2
  return(all(bool.vec))
}
D <- minimax(N=20,p=2,region="custom",const=ineqs)
plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2") #set up plot
polygon(c(0,2/3,1),c(0,2/3,0),col="gray") #design space
points(D,pch=16) #design points
mM <- mMdist(D,region="custom",const=ineqs)
mM$dist #minimax (fill) distance
lines(rbind(mM$far.pt,mM$far.despt),col="red",lty=2,lwd=2) #plot farthest point

#30-point minimax design on the hypercube [0,1]^6
D <- minimax(N=30,p=6)
mMdist(D)$dist #minimax (fill) distance

}
}
