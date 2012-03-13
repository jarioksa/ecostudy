\name{rmprey}
\alias{rmprey}
\alias{plot.rmprey}
\alias{traj.rmprey}
\alias{lines.rmprey}

\title{
  Rosenzweig-MacArthur Model of Predator-Prey Dynamics
}

\description{ The function sets up a Rosenzweig-MacArthur model of
  predator-prey relationship.  The prey has density-dependent control,
  and predator has type II functional response.  
}

\usage{
rmprey(b, K, w, D, e, s)
\method{plot}{rmprey}(x, xlim, ylim, arrows = 5, ...)
\method{lines}{rmprey}(x, N, P, time = 100, step = 0.2, ...)
\method{traj}{rmprey}(x, N, P, time = 100, step = 1, ...)
}

\arguments{
  \item{b, K, w, D, e, s}{Model parameters: see Details.}
  \item{x}{An \code{rmprey} result object.}
  \item{xlim, ylim}{Axis limits.}
  \item{arrows}{Number of trend arrows displayed in the graphs.}
  \item{N, P}{Initial prey and predator densities.}
  \item{time, step}{Total time span and time step for which
          trajectories are found.}
  \item{\dots}{Other parameters passed to graphical functions.}
}

\details{ 

  The prey \eqn{N} and predator \eqn{P} population changes are given
  by equations:
  \deqn{\frac{dN}{dt} = bN \left(1 - \frac{1}{K} N\right) - w \frac{N}{D+N}P}{dN/dt = b*N(1-N/K) - w*(N/(D+N))*P}
  \deqn{\frac{dP}{dt} = ew \frac{N}{D+N}P - sP}{dP/dt = e*w*(N/(D+N))*P - s*P}
  where \eqn{b} is the prey population growth rate, \eqn{K} is the
  prey carrying capacity, \eqn{w} is the attack rate, \eqn{D} is the
  handling time in in type II functional response, \eqn{e} is the
  production efficiency of the predator, and \eqn{s} is the predator
  mortality rate. Except for parameters \eqn{K} and \eqn{D} these are
  the same parameters as in the Lotka-Volterra predator-prey model
  (\code{\link{lotkaprey}}).  

  Function \code{rmprey} set ups the model, and \code{plot} draws the
  phase diagram.
}

\value{
  
  Function \code{rmprey} returns an object of class \code{"rmprey"}
  with following items:
  \item{preyfun}{A quadratic function giving the prey isocline.}

  \item{invpreyfun}{The inverse of \code{preyfun} that gives the (two)
  roots for given predator density, or \code{NA} if the roots are not
  real.}

  \item{preyopt, preytop, preylimits}{Descriptive statistics for the
  prey isocline: the location of the peak, the height of the peak, and
  the limits of the prey isocline. The limits are the real roots of
  \code{rmprey} or zero if the roots are negative.}

  \item{prediso}{The predator isocline.}

  \item{b, K, w, s, e, D}{Model parameters.}

  \item{call}{The function call.}
}

\author{ Jari Oksanen } 

\seealso{
\code{\link{lotkaprey}}.
}

\examples{
mod <- rmprey(0.8, 1000, 5, 400, 0.07, 0.2)
plot(mod)
lines(mod, 1000, 100)
plot(traj(mod, 1000, 100))
}

\keyword{ models }

