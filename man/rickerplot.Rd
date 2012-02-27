\name{rickerplot}
\alias{rickerplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rickerplot(N, R, K, b)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{N}{ ~~Describe \code{N} here~~ }
  \item{R}{ ~~Describe \code{R} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{b}{ ~~Describe \code{b} here~~ }
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
function (N, R, K, b) 
{
    i <- 1
    a <- ((R - 1)^(1/b))/K
    over <- 0
    tmp <- N
    repeat {
        N <- (N * R)/(1 + (a * N)^b)
        tmp <- c(tmp, N)
        i <- i + 1
        if (N > 0.99 * K) 
            over <- over + 1
        if (over > 10) 
            break
    }
    steps <- i
    plot(tmp, type = "l", xlab = "t", ylab = "N")
    abline(h = K, col = "red")
    x <- seq(from = 0, to = max(tmp) + 0.5, len = 100)
    plot(x, (x * R)/(1 + (a * x)^b), type = "l", xlab = expression(N[t]), 
        ylab = expression(N[t + 1]))
    lines(x, x, type = "l", lty = 2)
    for (i in 1:(steps - 1)) {
        now <- tmp[i]
        nex <- tmp[i + 1]
        segments(now, now, now, nex, col = "red")
        segments(now, nex, nex, nex, col = "red")
        i <- i + 1
    }
  }
}
\keyword{ ~keyword }%-- one or more ...
