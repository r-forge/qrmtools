\name{GEV_shape_plot}
\alias{GEV_shape_plot}
\title{Fitted GEV Shape as a Function of the Threshold}
\description{
  Fit GEVs to block maxima and plot the fitted GPD shape as
  a function of the block size.
}
\usage{
GEV_shape_plot(x, blocksize = tail(pretty(seq_len(length(x)/20), n = 64), -1),
               estimate.cov = TRUE, conf.level = 0.95,
               lines.args = list(lty = 2), xlab = "Block size",  ylab = NULL,
               xlab2 = "Number of blocks", ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data.}
  \item{blocksize}{\code{\link{numeric}} vector of block sizes for which
    to fit a GEV to the block maxima.}
  \item{estimate.cov}{\code{\link{logical}} indicating whether
    confidence intervals are to be computed.}
  \item{conf.level}{confidence level of the confidence intervals if
    \code{estimate.cov}.}
  \item{lines.args}{\code{\link{list}} of arguments passed to
    the underlying \code{\link{lines}()} for drawing the confidence intervals.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label (if \code{NULL}, a default is used).}
  \item{xlab2}{label of the secondary x-axis.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{
  \code{\link{invisible}()}
}
\details{
  Such plots can be used in the block maxima method for determining
  the optimal block size (as the smallest after which the plot
  is (roughly) stable).
}
\author{Marius Hofert}
\examples{
set.seed(271)
X <- rPar(1e5, shape = 4)
GEV_shape_plot(X)
abline(h = 1/4, lty = 3) # theoretical xi = 1/shape for Pareto
}
\keyword{parametric}
\keyword{hplot}