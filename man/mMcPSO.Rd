\name{mMcPSO}
\alias{mMcPSO}
\docType{package}
\title{
  Compute minimax designs using clustering
}
\description{
  \code{mMcPSO} is the main function for generating minimax designs on a desired design space \eqn{X} using a hybrid clustering - particle swarm optimization (PSO) algorithm. Subfunctions for \code{mMcPSO} are written in C++ for speed. Users have the flexibility of adjusting a variety of algorithmic parameters, including particle swarm optimization (PSO) settings, termination conditions, number of approximating points, etc. A formal exposition of this algorithm can be found in Mak and Joseph (2016).
}
\usage{
mMcPSO(N,p,q=10,
          pso=list(w=0.72,c1=1.49,c2=1.49),
          point_num=1e5,eval_num=10*point_num,point=NA,eval_pts=NA,
          bd = c(0,1),
          part_num=c(pso=10,pp=10),
          it_max=c(pso=200,pp=150,inn=1e4),
          it_lim=c(pso=25,pp=50),
          it_tol=c(pso=1e-4,pp=1e-4,inn=sqrt(p)*1e-4),
          region="uh",
          regionby=ifelse(p>2,1e-3,-1),
          jit=0.1/sqrt(N),
          pp_flag=F)
}
\arguments{
  \item{N}{Number of design points desired.}
  \item{p}{Dimension of design desired.}
  \item{q}{The approximation constant used to estimate the minimax criterion; refer to paper for details. Larger values of q give a better approximation, but may cause numerical instability.}
  \item{pso}{PSO settings for particle momentum (w), local-best velocity (c1) and global-best velocity (c2).}
  \item{point_num}{Number of points used to estimate the design space \eqn{X} for minimax clustering.}
  \item{eval_num}{Number of points used to estimate the design space \eqn{X} for post-processing.}
  \item{point,eval_pts}{User-specified representative points for clustering and post-processing. If NA, the algorithm generates these points using low-discrepancy sequences.}
  \item{bd}{Lower and upper bounds for the design space \eqn{X}.}
  \item{part_num}{Number of PSO particles for minimax clustering and post-processing.}
  \item{it_max}{Maximum iterations of minimax clustering, post-processing and the inner-loop for computing \eqn{C_q}-centers.}
  \item{it_lim,it_tol}{Algorithm terminates if the global-best objective does not improve by at least \code{it_tol} after \code{it_lim} iterations.}
  \item{region}{Option for non-hypercube design regions: Current choices include the unit hypercube "uh", the unit simplex "simp", and the unit ball "ball"}
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
plot(desuh,xlim=c(0,1),ylim=c(0,1))

#Generate and plot a minimax design with 20 points on the unit simplex A_2
# ... (CtoA provides the mapping from [0,1]^2 to A_2, see ?CtoA)
dessimp <- mMcPSO(N=20,p=2,region="simp")
plot(dessimp,xlim=c(0,1),ylim=c(0,1))
abline(0,1)

#Generate and plot a minimax design with 20 points on the unit ball B_2
# ... (CtoB2 provides the mapping from [0,1]^2 to B_2, see ?CtoB2)
library(plotrix)
desball <- mMcPSO(N=20,p=2,region="ball")
plot(desball,xlim=c(-1,1),ylim=c(-1,1))
draw.circle(0,0,1) #design boundaries

#Generate and plot a minimax design with 20 points on the unit hypercube [0,1]^4
desuh <- mMcPSO(N=20,p=4)
pairs(desuh,xlim=c(0,1),ylim=c(0,1))
}
}
