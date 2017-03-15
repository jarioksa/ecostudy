\encoding{UTF-8}
\name{ricker}
\alias{ricker}
\alias{plot.ricker}
\alias{traj.ricker}
 
\title{Ricker's Discrete Population Growth with Time-lag}

\description{

  Function finds a deterministic discrete population growth model
  (Ricker model) and draws population size against time and a phase
  diagram of population size and population increase side by side.

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
  \deqn{N_{t+1} = R N_t / (1 + (aN_t)^b) \,}{N[t+1] = r*N[t]/(1+(a*N[t])^b),}
  where
  \deqn{a = ((R-1)^{1/b})/K.}

  The current formulation follows Begon et al. (1996, p. 240--242).
  Krebs (2009, p. 153) uses different parametrization and discusses the
  model under \dQuote{time lag models}.

}

\value{
  Function returns the input model parameters, the calculated new
  parameter \eqn{a} and population sizes \eqn{N} at each step, and the
  function \code{call}.
}

\references{

  Begon, M., Harper, J.L. & Townsend, C.R. (1996) \emph{Ecology:
  Individuals, Populations and Communities.} Blackwell. 3 ed., 1068 p.

  Krebs, C. J. (2009) \emph{Ecology.} Benjamin Cummings. 6 ed., 655 p.

}

\author{Jari Oksanen}

\examples{
par(mfrow=c(2,2))
mod1 <- ricker(10, 2, 50, 1)
mod4 <- ricker(10, 2, 50, 4)
mod1
mod4
plot(mod1)
plot(mod4)
par(mfrow=c(1,1))
}
\keyword{ models }
