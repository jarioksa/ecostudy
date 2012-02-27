\name{poplogist}
\alias{poplogist}

\title{
  Logistic Population Growth
}
\description{
  Function returns the expected value at a given time for the logistic
  populatoin growth.
}
\usage{
poplogist(t, r, N0, K)
}

\arguments{
  \item{t}{
    Time.
}
  \item{r}{
    Instantaneous population growth rate.
}
  \item{N0}{
    Population size at time \eqn{t = 0}.
}
  \item{K}{
    Carrying capacity.
}
}
\details{
  \deqn{N = \frac{K}{1 + \frac{K-N_0}{N_0}\exp(-rt)}}{N = K/(1 + (K-N0)/N0 * exp(-r * t))}
}
\value{
  Function returns a vector of populations sizes corresponding to input
  vector of \code{t}.
}

\author{
  Jari Oksanen
}

\examples{
x <- seq(0,100)
plot(x, poplogist(x, 0.1, 10, 100), xlab = "Time (t)",
  ylab="Population size (N)", type = "l")
abline(h=100, lty=2)
}

\keyword{ models }
