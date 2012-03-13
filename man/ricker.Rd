\encoding{UTF-8}
\name{ricker}
\alias{ricker}
\alias{plot.ricker}
\alias{traj.ricker}
\alias{print.ricker}
 
\title{Ricker's Discrete Population Growth with Time-lag}

\description{

  Function finds a deterministic discrete population growth
  model (Ricker model) and draws side by side population size against
  time and a phase diagram of population size and population increase. 

}

\usage{
ricker(N, R, K, b)
\method{plot}{ricker}(x, which = c(1, 2), ...)
\method{traj}{ricker}(x, ...)
}

\arguments{
  \item{N}{Initial population size. }
  \item{R}{Population growth rate. }
  \item{K}{Carrying capacity. }
  \item{b}{Time lag parameter. }
  \item{x}{\code{ricker} result object.}
  \item{which}{Kind of graph displayed: \code{which = 1} plots
    population size against time using \code{\link{plot.traj}}, and
    \code{which = 2} displays the phase diagram with population
    trajectory.}
  \item{\dots}{Other parameters passed to graphical functions.}
}

\details{

  The population size \eqn{N} at time \eqn{t+1} is given as
  \deqn{N_{t+1} = R N_t / (1 + (aN)^b) \,}{N[t+1] = r*N[t]/(1+(a*N)^b),}
  where
  \deqn{a = ((R-1)^{1/b})/K.}
}

\value{
  Function returns the input model parameters, the calculated new
  parameter \eqn{a} and population sizes \eqn{N} at each step, and the
  function \code{call}.
}

\examples{
par(mfrow=c(2,2))
mod1 <- ricker(10, 2, 50, 1)
mod4 <- ricker(10, 2, 50, 4)
plot(mod1)
plot(mod4)
par(mfrow=c(1,1))
}
\keyword{ models }
