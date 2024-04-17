\name{gsummary}
\alias{gsummary}
\alias{gsummary.default}
\alias{gsummary.data.frame}
\alias{gsummary.lm}
\alias{gsummary.glm}

\title{The \code{gsummary()} function}

\description{The \code{gsummary()} function provides an alternative to the \code{\link{summary}()} function by returning different information. The prefix \code{g} is a reminder of who to blame if things do not work properly. }

\usage{
## generic:
gsummary(object, ...)
\method{gsummary}{default}(object, ...)
\method{gsummary}{data.frame}(object, ...)
\method{gsummary}{lm}(object, vcov.type = c("ordinary", "robust", "hac"), confint.level = 0.95, ...)
\method{gsummary}{glm}(object, confint.level = 0.95, ...)
}

\arguments{
  \item{object}{an object of suitable class, for example \code{\link{data.frame}}, \code{\link{lm}} or \code{\link{glm}}. }
  \item{vcov.type}{a character string that determines the variance-vcovariance estimator. If \code{"ordinary"} (default), then the ordinary estimator is used (\code{\link{vcov.lm}()}). If \code{"robust"}, then the heteroscedasticity robust estimator of White (1980) (\code{\link{vcovHC}()} with \code{type = "HC"}) is used. If \code{"hac"}, then the heteroscedasticity and autocorrelation robust estimator of Newey and West (1987) (\code{\link{NeweyWest}()}) is used. }
  \item{confint.level}{a number between 0 and 1 (the default is \code{0.95}), or \code{NULL}. If a number, then confidence intervals are printed (the default is 95 percent). If \code{NULL}, then confidence intervals are not printed. }
  \item{\dots}{additional arguments}
}

\value{No value is returned, the function only prints. The content of the print depends on the class of its main argument \code{object}. }

\references{
Halbert White (1980): 'A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity', Econometrica 48, pp. 817-838.
Whitney K. Newey and Kenned D. West (1987): 'A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix', Econometrica 55, pp. 703-708.
}

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{
\code{\link{summary}()}
}

\examples{
##simulate some data:
set.seed(123)
y <- rnorm(20); x <- rnorm(20); z <- rnorm(20)

##illustrate gsummary.data.frame():
mydataframe <- as.data.frame(cbind(y,x,z))
gsummary(mydataframe)

##illustrate gsummary.lm():
mymodel <- lm(y ~ x + z)
gsummary(mymodel)
gsummary(mymodel, vcov.type="robust")
gsummary(mymodel, vcov.type="hac")
gsummary(mymodel, confint.level=0.90)
gsummary(mymodel, confint.level=0.99)
gsummary(mymodel, confint.level=NULL)

##illustrate gsummary.glm():
y <- as.numeric( y > 0 )
mymodel <- glm(y ~ x + z, family=binomial)
gsummary(mymodel)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
