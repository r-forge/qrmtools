\name{edf_plot}
\alias{edf_plot}
\title{Plot Empirical Distribution Functions}
\description{
  Plotting empirical distribution functions.
}
\usage{
edf_plot(x, do.points = nrow(x) <= 100, log = "", xlim = NULL, ylim = NULL,
         col = NULL, main = "", xlab = "x", ylab = "Distribution function at x",
         ...)
}
\arguments{
  \item{x}{\code{\link{numeric}} vector of data (of which the empirical
    distribution function is to be plotted). Can also be a \code{\link{matrix}}
    in which case each column corresponds to one empirical distribution function.}
  \item{do.points}{\code{\link{logical}} indicating whether points are
    to be plotted; see \code{?plot.stepfun}.}
  \item{log}{\code{\link{character}} indicating whether
     a logarithmic x-axis is used (available are only \code{""} and \code{"x"}).}
  \item{xlim}{x-axis limits; defaults to the range of \code{x}.}
  \item{ylim}{y-axis limits; the lower bound defaults to 0, the upper to
    the maximal height of the empiricial distribution function(s).}
  \item{col}{color(s).}
  \item{main}{title.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{plot.stepfun}()}.}
}
\value{
  Nothing (plot by side-effect).
}
\author{Marius Hofert}
\examples{
x <- c(5, 2, 4, 2, 3, 2, 2, 2, 1, 2) # example data
edf_plot(x) # the default
edf_plot(x, verticals = FALSE) # the 'mathematical' version
edf_plot(x, do.points = FALSE) # good for many sample points
edf_plot(x, log = "x") # logarithmic; cannot show flat part before first jump
edf_plot(cbind(x, x+2), col = c("black", "royalblue3")) # plots both empirical distributions
edf_plot(cbind(x, x+2), log = "x", col = c("black", "royalblue3"))
}
\keyword{hplot}