\name{mMcPSO}
\alias{mMcPSO}
\docType{package}
\title{
  Compute minimax designs using clustering
}
\description{
  \code{mMcPSO} is the main function for generating minimax designs on a desired design space \eqn{X} using a hybrid clustering - particle swarm optimization (PSO) algorithm. Subfunctions for \code{mMcPSO} are written in C++ for speed. Users have the flexibility of adjusting a variety of algorithmic parameters, including particle swarm optimization (PSO) settings, termination conditions, number of approximating points, etc. A formal exposition of this algorithm can be found in Mak and Joseph (2017).
}
\usage{
mMcPSO(N,p,q=10,region="uh",
        pso=list(w=0.72,c1=1.49,c2=1.49),
        part_num_pso=10,part_num_pp=5,
        point_num=1e5,eval_num=10*point_num,point=NA,eval_pts=NA,
        it_max_pso=200,it_max_pp=ifelse(region=="simp",0,50),it_max_inn=1e4,
        it_lim_pso=25,it_lim_pp=10,
        it_tol_pso=1e-4,it_tol_pp=1e-4,it_tol_inn=1e-4,
        regionby=ifelse(p>2,1e-3,-1),
        jit=ifelse(region=="simp",0,0.1/sqrt(N)),
        pp_flag=F)
}
\arguments{
  \item{N}{Number of design points desired.}
  \item{p}{Dimension of design desired.}
  \item{q}{The approximation constant used to estimate the minimax criterion; refer to paper for details. Larger values of q give a better approximation, but may cause numerical instability.}
  \item{region}{Option for non-hypercube design regions: Current choices include the unit hypercube "uh", the unit simplex "simp", and the unit ball "ball"}
  \item{pso}{PSO settings for particle momentum (w), local-best velocity (c1) and global-best velocity (c2).}
  \item{part_num_pso,part_num_pp}{Number of PSO particles for minimax clustering and post-processing.}
  \item{point_num}{Number of points used to estimate the design space \eqn{X} for minimax clustering.}
  \item{eval_num}{Number of points used to estimate the design space \eqn{X} for post-processing.}
  \item{point,eval_pts}{User-specified representative points for clustering and post-processing. If NA, the algorithm generates these points using low-discrepancy sequences.}
  \item{it_max_pso,it_max_pp,it_max_inn}{Maximum iterations of minimax clustering, post-processing and the inner-loop for computing \eqn{C_q}-centers.}
  \item{it_lim_pso,it_lim_pp,it_tol_pso,it_tol_pp,it_tol_inn}{Algorithm terminates if the global-best objective does not improve by at least \code{it_tol} after \code{it_lim} iterations.}
  \item{regionby}{Specifies step-size for approximating integrals in non-hypercube transformations.}
  \item{jit}{Jitter radius for post-processing.}
  \item{pp_flag}{Redundant; still in development.}
}
\value{
  An \code{N}-by-\code{p} matrix representing the minimax design.
}
\examples{

\dontrun{
#Generate and plot a minimax design with 20 points on the unit hypercube [0,1]^2
desuh <- mMcPSO(N=20,p=2)
plot(desuh,xlim=c(0,1),ylim=c(0,1),pch=16)

#Generate and plot a minimax design with 20 points on the unit simplex A_2
# ... (CtoA provides the mapping from [0,1]^2 to A_2, see ?CtoA)
dessimp <- mMcPSO(N=20,p=2,region="simp")
plot(dessimp,xlim=c(0,1),ylim=c(0,1),pch=16)
abline(0,1)

#Generate and plot a minimax design with 20 points on the unit ball B_2
# ... (CtoB2 provides the mapping from [0,1]^2 to B_2, see ?CtoB2)
library(plotrix)
desball <- mMcPSO(N=20,p=2,region="ball")
plot(desball,xlim=c(-1,1),ylim=c(-1,1),pch=16)
draw.circle(0,0,1) #design boundaries

#Generate and plot a minimax design with 20 points on the unit hypercube [0,1]^4
desuh <- mMcPSO(N=20,p=4)
pairs(desuh,xlim=c(0,1),ylim=c(0,1),pch=16)
}
}
