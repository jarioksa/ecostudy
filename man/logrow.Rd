\name{logrow}
\alias{logrow}

\title{Random Logistic Population Growth}
\description{
  Random population growth with logistic expectation.
}
\usage{
logrow(N0, b, d, K, bcontr=0.5, out=0)
}

\arguments{
  \item{N0}{Initial population size }
  \item{b}{Birth rate. }
  \item{d}{Death rate. }
  \item{K}{Carrying capacity. }
  \item{bcontr}{The proportion of decrease of growth rate due to
    decrease in the birth rate. }
  \item{out}{Desire output? }
}
\details{
  The stochastic population growth rate is based on logistic expectation
  following Pielou's idea.
}
\value{
  Output can be suppressed and then only a graph is drawn. The possible
  output is
  \item{time }{time.}
  \item{simN }{simulated population size.}
}

\author{Jari Oksanen }

\examples{
par(mfrow=c(2,2))
replicate(4, logrow(5, 1.2, 0.8, 50))
par(mfrow=c(1,1))
}
\keyword{ ~keyword }%-- one or more ...
