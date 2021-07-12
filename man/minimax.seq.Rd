\name{minimax.seq}
\alias{minimax.seq}
\docType{package}
\title{
  Compute sequential minimax designs using clustering on constrained design regions
}
\description{
  \code{minimax.seq} computes sequential minimax designs using the minimax clustering algorithm in Mak and Joseph (2018). Current design region options include the unit hypercube ("hypercube"), the unit simplex ("simplex"), the unit ball ("ball"), a user-defined active subspace ("subspace"), and user-defined constraints on the unit hypercube ("custom").
}
\usage{
minimax.seq(D,Nseq,q=10,region="hypercube",const=NA,subsp=NA,wt.subsp=FALSE,
          params_pso=list(w=0.72,c1=1.49,c2=1.49),
          npart=5,nclust=1e5,neval=nclust,
          itmax_pso=50,itmax_pp=100,itmax_inn=1e4,jit=0.1/sqrt(Nseq+nrow(D)))
}
\arguments{
  \item{D}{Initial design.}
  \item{Nseq}{Number of sequential design points.}
  \item{q}{Power parameter for approximating the minimax criterion (see paper for details). Larger values of q give a better approximation, but may cause numerical instability.}
  \item{region}{String for desired design region. Current options include "hypercube", "simplex", "ball", "subspace", and "custom".}
  \item{const}{Function for desired constraints (inequalities) on design space. See examples for implementation.}
  \item{subsp}{A \code{p}' x \code{p} matrix for the desired active subspace, where \code{p}' and \code{p} are the dimensions of the full space and active subspace. Only used when \code{region} is "subspace".}
  \item{wt.subsp}{TRUE if a weighted design is desired over the active subspace, FALSE if an unweighted (uniform) is desired. Only used when \code{region} is "subspace".}
  \item{params_pso}{Particle swarm optimization parameters (particle momentum (w), local-best velocity (c1) and global-best velocity (c2)).}
  \item{npart}{Number of particles for particle swarm optimization.}
  \item{nclust,neval}{Number of sample points for minimax clustering and post-processing.}
  \item{itmax_pso,itmax_pp,itmax_inn}{Maximum number of iterations for minimax clustering, post-processing and inner optimization.}
  \item{jit}{Jitter radius for post-processing.}
}
\value{
  An \code{Nseq}-by-\code{p} matrix for the minimax design.
}
\examples{

\dontrun{
# 20+20-point sequential minimax design on the hypercube [0,1]^2
Nini <- 20
Nseq <- 20
Dini <- matrix(runif(Nini*2),ncol=2) # random initial design
Dseq <- minimax.seq(Dini,Nseq) # sequential minimax design

plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2") #set up plot
polygon(c(0,0,1,1),c(0,1,1,0),col="gray") #design space
points(Dini,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2",pch=16) #original design
points(Dseq,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2",pch=16,col="blue") #sequential design

# 20+20-point sequential minimax design on active subspace
set.seed(1)
library(rstiefel)
library(geometry)
p.full <- 10 #full dimension
p.as <- 2 #active subspace dimension
Nini <- 20
Nseq <- 20
Dini <- matrix(runif(Nini*p.full),ncol=p.full) # random initial design
A <- rustiefel(p.full,p.as) # random active subspace
Dini <- t(t(A)\%*\%t(Dini)) # initial design projected on active subspace
# sequential minimax design
Dseq <- minimax.seq(Dini,Nseq=20,region="subspace",subsp=A,wt.subsp=TRUE)

# Compute active subspace polygon
Nplot <- 1e5
plot.pts <- matrix(runif(Nplot*p.full),ncol=p.full)
plot.pts <- t(t(A)\%*\%t(plot.pts))
hull <- convhulln(plot.pts)
for (i in 2:nrow(hull)){
  newidx <- setdiff(c(which(hull[,2]==hull[i-1,2]),which(hull[,1]==hull[i-1,2])),i-1)
  if (hull[newidx,1]==hull[i-1,2]){
    tmp <- hull[newidx,]
    hull[newidx,] <- hull[i,]
    hull[i,] <- tmp
  }else{
    tmp <- rev(hull[newidx,])
    hull[newidx,] <- hull[i,]
    hull[i,] <- tmp
  }
}

# Visualize
plot(NULL,xlim=c(-2,2),ylim=c(-2,2),xlab="x1",ylab="x2") #set up plot
polygon(plot.pts[hull[,1],],lwd=2,col="gray") #plot active subspace
points(Dini,xlab="x1",ylab="x2",pch=16) #original design
points(Dseq,xlab="x1",ylab="x2",pch=16,col="blue") #sequential design
}
}
