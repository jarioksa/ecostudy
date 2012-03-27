\name{ranlogist}
\alias{ranlogist}
\alias{plot.ranlogist}
\alias{lines.ranlogist}
\alias{traj.ranlogist}
\alias{print.ranlogist}

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

  The stochastic population growth rate is based on logistic
  expectation following Pielou (1969). The instantaneous growth rate
  is the difference of birth and death rates \eqn{r = b - d}.  The
  birth and death rates are dependent on current population size and
  at carrying capacity \eqn{K} they are equal (\eqn{r = b}). The sum
  of death and birth rates at given population size gives the current
  event rate which is used to simulate the time to the next event from
  exponential distribution.  The next event is either birth that
  increases population size by one or a death that reduces population
  size by one individual. Consequently, the population size is
  discrete but time is continuous, unlike in many other stochastic
  models which treat time as discrete, but population size as
  continuous.  Further, Pielou's model has real temporal
  autocorrelation of population size, unlike many other models which
  vary expected population growth with no autocorrelation.

  Krebs (2009, p. 155) discusses Pielou (1969) briefly, but ignores
  the advanced features of her model.

  Function \code{ranlogist} perform the simulation, \code{traj}
  extracts the simulation results, and \code{plot} and \code{lines}
  can be used for graphics.  The \code{plot} is based on
  \code{\link{plot.traj}}.

}

\value{
  Function \code{ranlogist} returns the model parameters, \eqn{r},
  targeted end time of simulation (\code{timend}) and the simulation
  results in vectors \code{time} and \code{N}, and the function
  \code{call}. 
}

\references{
  
  Krebs, C. J. (2009) \emph{Ecology.} Benjaming Cummings. 6 ed., 655 pp.

  Pielou, E.C. (1969) \emph{An Introduction to Mathematical Ecology.}
  Wiley-Interscience, 286 pp.
 
}

\author{Jari Oksanen }

\examples{
par(mfrow=c(2,2))
replicate(4, plot(ranlogist(5, 1.2, 0.8, 50)))
par(mfrow=c(1,1))
}
\keyword{ models }
