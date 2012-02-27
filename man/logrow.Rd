\name{logrow}
\alias{logrow}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logrow(N0, b, d, K, bcontr=0.5, out=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{N0}{ ~~Describe \code{N0} here~~ }
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{d}{ ~~Describe \code{d} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{bcontr}{ ~~Describe \code{bcontr} here~~ }
  \item{out}{ ~~Describe \code{out} here~~ }
}
\details{
  ~~ If necessary, more details than the __description__  above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of `comp1'}
  \item{comp2 }{Description of `comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	     or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (N0, b, d, K, bcontr = 0.5, out = 0) 
{
    c <- (K - N0)/N0
    r <- b - d
    t <- 0
    N <- N0
    while (N[t + 1] < 0.998 * K) {
        t <- t + 1
        tmp <- K/(1 + c * exp(-r * t))
        N <- c(N, tmp)
    }
    tend <- length(N)
    t <- seq(0, length(N) - 1)
    plot(t, N, type = "l", col = "red", ylim = c(0, 1.2 * K))
    abline(h = K)
    time <- 0
    now <- 0
    simN <- N0
    size <- N0
    while (now < tend) {
        adj <- (b - d) * (size/K)
        bN <- b - bcontr * adj
        dN <- d + (1 - bcontr) * adj
        events <- (bN + dN) * size
        step <- rexp(1, events)
        now <- now + step
        if (runif(1) <= bN/(bN + dN)) {
            size <- size + 1
        }
        else {
            size <- size - 1
        }
        simN <- c(simN, size)
        time <- c(time, now)
        if (size == 0) 
            break
    }
    lines(time, simN, col = "blue")
    final <- cbind(time, simN)
    if (out != 0) 
        return(final)
  }
}
\keyword{ ~keyword }%-- one or more ...
