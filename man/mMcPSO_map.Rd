\name{mMcPSO_map}
\alias{mMcPSO_map}
\docType{package}
\title{
  Compute minimax designs using clustering on a user-provided map (provided as an image file).
}
\description{
  \code{mMcPSO_map} generates minimax designs on a user-provided map (provided as an image file) using a hybrid clustering - particle swarm optimization (PSO) algorithm. This function uses the minimax clustering routine \code{mMcPSO} internally as a workhorse. A formal exposition of this algorithm can be found in Mak and Joseph (2017).
}
\usage{
mMcPSO_map(N,img,p=2,q=10,
            pso=list(w=0.72,c1=1.49,c2=1.49),
            part_num_pso=10,part_num_pp=5,
            point_num=1e5,eval_num=10*point_num,point=NA,eval_pts=NA,
            it_max_pso=200,it_max_pp=50,it_max_inn=1e4,
            it_lim_pso=25,it_lim_pp=10,
            it_tol_pso=1e-4,it_tol_pp=1e-4,it_tol_inn=1e-4,
            jit=0.1/sqrt(N))
}
\arguments{
  \item{N}{Number of design points desired.}
  \item{img}{A 0-1 matrix corresponding to the image file of the map.}
  \item{p}{Dimension of design desired (default = 2).}
  \item{q}{The approximation constant used to estimate the minimax criterion; refer to paper for details. Larger values of q give a better approximation, but may cause numerical instability.}
  \item{pso}{PSO settings for particle momentum (w), local-best velocity (c1) and global-best velocity (c2).}
  \item{part_num_pso,part_num_pp}{Number of PSO particles for minimax clustering and post-processing.}
  \item{point_num}{Number of points used to estimate the design space \eqn{X} for minimax clustering.}
  \item{eval_num}{Number of points used to estimate the design space \eqn{X} for post-processing.}
  \item{point,eval_pts}{User-specified representative points for clustering and post-processing. If NA, the algorithm generates these points using low-discrepancy sequences.}
  \item{it_max_pso,it_max_pp,it_max_inn}{Maximum iterations of minimax clustering, post-processing and the inner-loop for computing \eqn{C_q}-centers.}
  \item{it_lim_pso,it_lim_pp,it_tol_pso,it_tol_pp,it_tol_inn}{Algorithm terminates if the global-best objective does not improve by at least \code{it_tol} after \code{it_lim} iterations.}
  \item{jit}{Jitter radius for post-processing.}
}
\value{
  An \code{N}-by-\code{p} matrix representing the minimax design.
}
\examples{
  \dontrun{
    #Generate and plot a minimax design with 20 points on the map of Georgia
    library(jpeg)
    n <- 25
    img <- readJPEG(system.file("img", "gamap.jpg", package="minimaxdesign"))[,,1]
    image(t(img)[,nrow(img):1],col=gray.colors(12,start=0.6),main="Georgia")
    img <- t(img)[,nrow(img):1] #Invert image due to reading distortion
    des <- mMcPSO_map(n,img)
    points(des,pch=16)
  }
}
