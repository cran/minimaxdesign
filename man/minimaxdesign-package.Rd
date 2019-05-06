\name{minimaxdesign-package}
\alias{minimaxdesign-package}
\alias{minimaxdesign}
\docType{package}
\title{An R package for computing Minimax and Minimax Projection Designs}
\description{
The 'minimaxdesign' package provides functions for generating minimax designs and minimax projection designs.
}
\details{
\tabular{ll}{
Package: \tab minimaxdesign\cr
Type: \tab Package\cr
Version: \tab 0.1.4\cr
Date: \tab 2019-05-03\cr
License: \tab GPL (>= 2)\cr
}
Provides two main functions, \code{minimax()} and \code{miniMaxPro()}, for computing minimax     and minimax projection designs using the minimax clustering algorithm in Mak and
    Joseph (2018) <DOI:10.1080/10618600.2017.1302881>. Current design region options
    include the unit hypercube ("hypercube"), the unit simplex ("simplex"), the unit ball
    ("ball"), as well as user-defined constraints on the unit hypercube ("custom"). Minimax
    designs can also be computed on user-provided images using the function \code{minimax.map()}. Design quality can be assessed using the function \code{mMdist()}, which computes the minimax (fill) distance of a design.
}
\author{
Simon Mak

Maintainer: Simon Mak <smak6@gatech.edu>
}
\references{
Mak, S. and Joseph, V. R. (2018). Minimax and minimax projection designs using clustering. \emph{Journal of Computational and Graphical Statistics}, 27(1):166-178.
}

\keyword{ package }

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
  }

}
