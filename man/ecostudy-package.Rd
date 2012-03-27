\name{ecostudy-package}
\alias{ecostudy-package}
\alias{ecostudy}
\docType{package}
\title{
Ecology teaching package
}
\description{  
  
  Functions to aid understanding ecological models for Ecology 101
  level. The package is designed to support Krebs (2009), but because
  it covers the most basic models, it probably suits for any other
  similar textbook. The primary usage was to draw graphics for my
  basic ecology lectures, but many functions can be used for other
  teaching or studying purposes.

}
 
\details{

  The package contains basic ecological models for population dynamics
  of single species and two-species interactions.  The models
  include:
  \itemize{
  
  \item Population growth: logistic model (\code{\link{poplogist}}),
  stochastic logistic growth (\code{\link{ranlogist}}), Ricker model
  (\code{\link{ricker}}), Milne's stochastic ceiling model
  (\code{\link{ranmilne}}).

  \item Competition: Lotka-Volterra model (\code{\link{lotkacomp}}),
  Tilman's resource competition model (\code{\link{tilman}}).

  \item Predator--Prey relationship: Lotka-Volterra model in its basic
  form (\code{\link{lotkaprey}}), and the model with density dependent
  prey (\code{\link{lotkapreyK}}), Rosenzweig--MacArthur model
  (\code{\link{rmprey}}).  
  }

  The functions were originally written to draw graphics for my
  lecture slides, and therefore their user interface often is
  restricted (contributions are welcome). Most functions work so that
  the the basic function sets up the model, \code{plot} draws the zero
  isoclines and other static structures but no population
  trajectories. The population trajectories can be found with a new S3
  generic function \code{\link{traj}} in the package, and \code{lines}
  function normally adds the trajectories to the basic plot.

  Most functions are based on the \pkg{primer} package (Stevens 2009),
  and even other functions are modelled similarly. The models are set
  up as differential equations, and numerical integration is used to
  find the trajectories (function \code{\link[deSolve]{ode}} of the
  \pkg{deSolve} package).

}

\author{
Jari Oksanen <jari.oksanen@oulu.fi>

Maintainer: Jari Oksanen <jari.oksanen@oulu.fi>

}
\references{
Krebs, C. J. (2009) \emph{Ecology.} Benjamin Cummings. 6 ed., 655 p.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R.} Springer. 
}
\keyword{ package }
\seealso{
\code{\link[primer:primer-package]{primer}}.
}
\examples{
## A typical usage: Lotka-Volterra competition
mod <- lotkacomp(1.2, 1.1, 50, 50)
## The static model
mod
summary(mod)
plot(mod)
## Species trajectories
tr <- traj(mod, N1 = 20, N2 = 4, r1 = 0.3, r2 = 1)
plot(tr)
## Plot species trajectories
plot(mod)
lines(mod, N1 = 20, N2 = 4, r1 = 0.2, r2 = 1)
}
