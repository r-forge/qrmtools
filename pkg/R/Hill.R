### Hill estimator #############################################################

##' @title Hill Estimator
##' @param x vector of numeric data
##' @param k vector of length 2, determining the smallest and largest number
##'        of order statistics of x to compute the Hill estimator for (the
##'        smallest needs to be >= 2). If of length 1, k is expanded
##'        by length(x).
##' @param conf.level confidence level for the asymptotic CIs
##' @return (k[1]:k[2], 5)-matrix with columns:
##'         - k = the k's used
##'         - k.prob = empirical probabilities corresponding to the k's
##'         - tail.index = estimated tail indices for the different k's
##'         - CI.low = lower CI for the tail indices
##'         - CI.up  = upper CI for the tail indices
##' @author Marius Hofert
Hill_estimator <- function(x, k = c(10, length(x)), conf.level = 0.95)
{
    ## Basics
    if(!is.vector(x)) x <- as.vector(x) # convert time series to vector (otherwise sort() fails)
    n <- length(x)
    if(length(k) == 1) k <- c(k, n)
    stopifnot(is.numeric(x), x > 0, length(k) == 2, 2 <= k, diff(k) > 0,
              0 < conf.level, conf.level < 1)
    if(k[2] > n) # more meaningful error message
        stop("k[2] = ",k[2]," must be <= length(x) = ",n)

    ## Build ingredients
    lx <- sort(log(x), decreasing = TRUE) # sorted in decreasing order

    ## Hill estimator
    k.min <- k[1]
    k.max <- k[2]
    k.ran <- k.min:k.max
    k.prob <- 1 - (k.ran-1) / n # k.ran == 2 corresponds to 1 - 1/n, k.ran = n to 1/n
    one.to.k.max <- seq_len(k.max)
    lx.cmean.k.max <- cumsum(lx[one.to.k.max]) / one.to.k.max # compute all cumulative means from 1 to k.max
    tail.index <- 1 / (lx.cmean.k.max[k.ran] - lx[k.ran])

    ## This form is frequently found in the literature (although
    ## EKM (1997, p. 190, Eq. (4.12)) use the one with k ~> k-1)

    ## CI (see evir::hill or QRM::hillPlot)
    SE <- tail.index / sqrt(k.ran)
    q <- qnorm(1-(1-conf.level)/2)
    CI.low <- tail.index - q * SE
    CI.up  <- tail.index + q * SE

    ## Return
    cbind(k = k.ran, k.prob = k.prob, tail.index = tail.index,
          CI.low = CI.low, CI.up = CI.up)
}


### Hill plot ##################################################################

##' @title Hill Plot (Estimated Tail Index as a Function of the Number
##'        of Order Statistics)
##' @param x see ?Hill_estimator
##' @param k see ?Hill_estimator
##' @param conf.level see ?Hill_estimator
##' @param Hill.estimator object as returned by Hill_estimator()
##' @param log see ?plot
##' @param xlim see ?plot
##' @param ylim see ?plot
##' @param xlab see ?plot
##' @param ylab see ?plot
##' @param lines.args arguments passed to the underlying lines()
##' @param CI.col color of the confidence interval region
##' @param third.axis logical indicating whether the 3rd axis is
##'        plotted
##' @param tlab label of the third axis
##' @param ... additional arguments passed to the underlying plot()
##' @return Hill plot (by side-effect)
##' @author Marius Hofert
Hill_plot <- function(x, k = c(10, length(x)), conf.level = 0.95, Hill.estimator = NULL,
                      log = "x", xlim = NULL, ylim = NULL,
                      xlab = "Order statistics", ylab = "Tail index",
                      lines.args = list(),
                      CI.col = adjustcolor(1, alpha.f = 0.2),
                      third.axis = TRUE, tlab = "Empirical probability",
                      ...)
{
    ## Ingredients
    if(is.null(Hill.estimator))
        Hill.estimator <- Hill_estimator(x, k = k, conf.level = conf.level)
    k <- Hill.estimator[,"k"]
    k.prob <- Hill.estimator[,"k.prob"]
    tail.index <- Hill.estimator[,"tail.index"]
    CI.low <- Hill.estimator[,"CI.low"]
    CI.up  <- Hill.estimator[,"CI.up"]
    if(is.null(xlim)) xlim <- range(k)
    if(is.null(ylim))
        ylim <- if(is.na(CI.col)) range(tail.index) else range(tail.index, CI.low, CI.up)

    ## Hill plot with CI
    plot(NA, log = log, xlim = range(k), ylim = ylim, xlab = xlab, ylab = ylab, ...) # setup plot region
    polygon(x = c(k, rev(k)), y = c(CI.low, rev(CI.up)), col = CI.col, border = NA) # CI
    do.call(lines, args = c(list(x = k, y = tail.index), lines.args)) # Hill plot

    ## Third axis
    if(third.axis) {
        at <- axTicks(1) # as for x-axis
        labs. <- if(at[1] == 0) {
                     c(1, k.prob[at])
                 } else {
                     k.prob[at]
                 }
        labs <- format(signif(labs., 4))
        axis(3, at = at, labels = labs)
        mtext(tlab, side = 3, line = 3)
    }
}
