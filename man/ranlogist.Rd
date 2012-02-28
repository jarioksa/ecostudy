\name{ranlogist}
\alias{ranlogist}
\alias{plot.ranlogist}
\alias{lines.ranlogist}
\alias{traj.ranlogist}

\title{Random Logistic Population Growth}
\description{
  Random population growth with logistic expectation.
}
\usage{
ranlogist(N0, b, d, K, bcontr=0.5)
\method{plot}{ranlogist}(x, ...)
\method{lines}{ranlogist}(x, ...)
\method{traj}{ranlogist}(x, ...)
}

\arguments{
  \item{N0}{Initial population size }
  \item{b}{Birth rate. }
  \item{d}{Death rate. }
  \item{K}{Carrying capacity. }
  \item{bcontr}{The proportion of decrease of growth rate due to
    decrease in the birth rate. }
  \item{x}{\code{"ranlogist"} result object. }
  \item{\dots}{Other parameters passed to graphical functions.}
}
\details{
  The stochastic population growth rate is based on logistic expectation
  following Pielou's idea. Function \code{ranlogist} perform the
  simulation, \code{traj} extracts the simulation results, and
  \code{plot} and \code{lines} can be used for graphics.  The
  \code{plot} is based on \code{\link{plot.traj}}. 
}
\value{
  Function \code{ranlogist} returns the model parameters, \eqn{r},
  targeted end time of simulation (\code{timend}) and the simulation
  results in vectors \code{time} and \code{N}.
}

\author{Jari Oksanen }

\examples{
par(mfrow=c(2,2))
replicate(4, plot(ranlogist(5, 1.2, 0.8, 50)))
par(mfrow=c(1,1))
}
\keyword{ models }
