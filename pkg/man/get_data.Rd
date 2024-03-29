\name{get_data}
\alias{get_data}
\title{Tools for Getting and Working with Data}
\description{
  Download (and possibly) merge data from freely available databases.
}
\usage{
get_data(x, from = NULL, to = NULL,
         src = c("yahoo", "quandl", "oanda", "FRED", "google"),
         FUN = NULL, verbose = TRUE, warn = TRUE, ...)
}
\arguments{
  \item{x}{vector of ticker symbols (e.g. \code{"^GSPC"} if
    \code{src = "yahoo"} or \code{"EUR/USD"} if \code{src = "oanda"}).}
  \item{from}{start date as a \code{Date} object or character string
    (in international date format \code{"yyyy-mm-dd"}); if \code{NULL},
    the earliest date with available data is picked.}
  \item{to}{end date as a \code{Date} object or character string
    (in international date format \code{"yyyy-mm-dd"}); if \code{NULL},
    the last date with available data is picked.}
  \item{src}{character string specifying the data source
    (e.g. \code{"yahoo"} for stocks or \code{"oanda"} for FX data);
    see \code{\link{getSymbols}()} and \code{\link{Quandl}()}.}
  \item{FUN}{\code{\link{function}} to be applied to the data before
    being returned. This can be
    \describe{
      \item{the identity:}{if the data could not be retrieved (and is
	thus replaced by \code{\link{NA}});}
      \item{the given \code{FUN}:}{if \code{FUN} has been provided;}
      \item{a useful default:}{if \code{FUN = NULL}; the default uses
	the adjusted close price \code{\link{Ad}()} if
	\code{src = "yahoo"},
	the close price \code{\link{Cl}()} if
	\code{src = "google"} and the identity otherwise.}
    }}
  \item{verbose}{\code{\link{logical}} indicating whether progress
    monitoring should be done.}
  \item{warn}{\code{\link{logical}} indicating whether a warning
    is given showing the error message when fetching \code{x} fails.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{getSymbols}()} from \pkg{quantmod} or
    \code{\link{Quandl}()} from \pkg{Quandl} (if \code{src = "quandl"}).}
}
\value{
  \code{xts} object containing the data with column name(s)
  adjusted to be the ticker symbol (in case lengths match; otherwise the
  column names are not adjusted); \code{\link{NA}} if
  data is not available.
}
\details{
  \code{FUN} is typically one of \pkg{quantmod}'s \code{\link{Op}},
  \code{\link{Hi}}, \code{\link{Lo}}, \code{\link{Cl}}, \code{\link{Vo}},
  \code{\link{Ad}} or one of the combined functions \code{\link{OpCl}},
  \code{\link{ClCl}}, \code{\link{HiCl}}, \code{\link{LoCl}}, \code{\link{LoHi}},
  \code{\link{OpHi}}, \code{\link{OpLo}}, \code{\link{OpOp}}.
}
\author{Marius Hofert}
\examples{
\dontrun{
## Note: This needs a working internet connection
## Get stock and volatility data (for all available trading days)
dat <- get_data(c("^GSPC", "^VIX")) # note: this needs a working internet connection
## Plot them (Alternative: plot.xts() from xtsExtra)
library(zoo)
plot.zoo(dat, screens = 1, main = "", xlab = "Trading day", ylab = "Value")
}
}
\keyword{manip}