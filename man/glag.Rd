\name{glag}
\alias{glag}
\alias{glag.default}

\title{ Lag a vector or a matrix, with special treatment of \code{\link[zoo]{zoo}} and \code{\link[stats]{ts}} objects}

\description{
The \code{glag()} function is similar to the \code{\link[stats]{lag}()} function from the \pkg{stats} package, but \code{glag()} actually \emph{lags} (the default in \code{\link[stats]{lag}()} is to lead). The funtion \code{glag()} also enables padding (for example \code{NA}s or \code{0}s) of the lost entries. Contrary to the \code{\link[stats]{lag}()} function, however, the default in \code{glag()} is to pad (with \code{NA}s). The \code{glag()} is particularly suited for \code{\link[zoo]{zoo}} objects, since their indexing is retained. The prefix \code{g} is a reminder of who to blame if things do not work properly.
}

\usage{
## generic:
glag(x, ...)
\method{glag}{default}(x, k = 1, pad = TRUE, pad.value = NA, ...)
}

\arguments{
  \item{x}{a numeric vector or matrix, or objects of class \code{\link[zoo]{zoo}} or \code{\link[stats]{ts}}. }
  \item{k}{integer equal to the lag (the default is 1). Negative values (that is, 'leading') is not possible. }
  \item{pad}{logical. If \code{TRUE} (default), then the lost entries are padded with \code{pad.value}. If \code{FALSE}, then no padding is undertaken. }
  \item{pad.value}{the padding value. }
  \item{\dots}{additional arguments}
}

\value{A vector or matrix, or objects of class \code{\link[zoo]{zoo}} or \code{\link[stats]{ts}}, with the lagged values. }

\author{Genaro Sucarrat, \url{https://www.sucarrat.net/}}

\seealso{
\code{\link[stats]{lag}()}, \code{\link[zoo]{lag.zoo}()}
}

\examples{
##generate some data:
x <- rnorm(5)

##lag series with NAs on missing entries:
glag(x)

##lag series with no padding:
x <- rnorm(5)
glag(x, pad = FALSE)

##lag series and retain the original zoo-index ordering:
x <- as.zoo(rnorm(5))
glag(x)

##lag two periods:
glag(x, k = 2)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
