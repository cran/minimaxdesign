\name{minimax.map}
\alias{minimax.map}
\docType{package}
\title{
  Compute minimax designs using clustering on a user-provided image
}
\description{
  \code{minimax.map} computes minimax designs on a user-provided binary (0-1) image, using the minimax clustering algorithm in Mak and Joseph (2018).
}
\usage{
minimax.map(N,img,p=2,q=10,
            params_pso=list(w=0.72,c1=1.49,c2=1.49),
            npart=5,nclust=1e5,neval=nclust,
            itmax_pso=50,itmax_pp=100,itmax_inn=1e4,jit=0.1/sqrt(N))
}
\arguments{
  \item{N}{Number of design points.}
  \item{img}{A binary 0-1 matrix, with 1 indicating the desired design region.}
  \item{p}{Dimension of design region.}
  \item{q}{Power parameter for approximating the minimax criterion (see paper for details). Larger values of q give a better approximation, but may cause numerical instability.}
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
    library(jpeg)
    n <- 25
    img <- readJPEG(system.file("img", "gamap.jpg", package="minimaxdesign"))[,,1]
    image(t(img)[,nrow(img):1],col=gray.colors(12,start=0.6),main="Georgia")
    img <- t(img)[,nrow(img):1] #Invert image due to reading distortion
    des <- minimax.map(n,img)
    points(des,pch=16)
  }
}
