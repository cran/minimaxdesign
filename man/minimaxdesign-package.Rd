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
Version: \tab 1.0\cr
Date: \tab 2016-04-13\cr
License: \tab GPL (>= 2)\cr
}
The 'minimaxdesign' package provides two main functions: \code{mMcPSO()} and
    \code{miniMaxPro()}, which generates minimax designs and minimax projection designs using
    clustering and particle swarm optimization (PSO) techniques. These designs can be used
    in a variety of settings, e.g., as space-filling designs for computer experiments or
    sensor allocation designs. A detailed description of the two designs and the employed
    algorithms can be found in Mak and Joseph (2016).
}
\author{
Simon Mak

Maintainer: Simon Mak <smak6@gatech.edu>
}
\references{
Mak, S., & Joseph, V.R. (2016). Minimax designs using clustering. Journal of Computational and Graphical Statistics. In revision.
}

\keyword{ package }

\examples{

\dontrun{

#Generate and plot a minimax design with 7 points on the unit hypercube [0,1]^2
mMdes <- mMcPSO(N=7,p=2)

#Generate a miniMaxPro design of 20 points on the unit hypercube [0,1]^8
library(foreach)
mMPdes <- miniMaxPro(N=20,p=8)$miniMaxPro

}

}
