#-------------------------------------------
# Documentation of methods
#' @name methods-conjugate
#' @title S3 Methods for the Object of Conjugate Class
#' @description Methods functions
#' @param x An object of class \code{conjugate}.
#' @param object An object of class \code{conjugate}.
#' @param dist The distribution that is used, \code{posterior},
#'     \code{prior} or \code{likelihood} distribution.
#' @param ... unused.
#'
NULL

#-------------------------------------------
# Print
#' @rdname methods-conjugate
#' @export
print.conjugate <- function(x, ...) {
    dist <- gsub("^(.*)-.*$", "\\1", x$model)
    parts <- c("prior", "likelihood", "posterior")
    nc <- max(nchar(as.integer(unlist(x[parts]))))
    pars <- lapply(parts, function(m) {
        pars <- formatC(unlist(x[[m]]), digits = 2,
                        width = nc + 3, format = "f")
        pars <- paste(pars, collapse = ", ")
        paste0("(", pars, ")")
    })
    cat(sprintf("Conjugate distribution: %s\n", x$model), sep = "\n")
    cat(paste0("  Prior:      ", dist, pars[[1]]), sep = "\n")
    cat(paste0("  Likelihood: ", dist, pars[[2]]), sep = "\n")
    cat(paste0("  Posterior:  ", dist, pars[[3]]), sep = "\n")
    cat("", sep = "\n")
    invisible()
}

#-------------------------------------------
# Summary
#' @rdname methods-conjugate
#' @export
summary.conjugate <- function(object, ...) {
    dist <- c("posterior", "likelihood", "prior")
    probs <- c(0.025, 0.975)
    cname <- paste0("QTS", probs)
    conf_qts <- lapply(dist, function(d)
        confint_qts(object, dist = d))
    out <- do.call(rbind, object$summary[dist])
    out <- cbind(out, do.call(rbind, conf_qts))
    out <- matrix(unlist(out), nrow = 3)
    colnames(out) <- c("Mean", "Mode", "Variance", cname)
    rownames(out) <- dist
    attr(out, "class") <- "summary.conjugate"
    attr(out, "model") <- object$model
    return(out)
}

#-------------------------------------------
# Print summary
#' @rdname methods-conjugate
#' @export
#' @param digits the desired number of digits after the decimal point.
print.summary.conjugate <- function(x, digits = getOption("digits"),
                                    ...) {
    model <- attr(x, "model")
    nc <- max(nchar(as.integer(x)[!is.na(x)]))
    x <- formatC(x, digits = digits, width = nc + digits + 1,
                 format = "f")
    conf_qts <- apply(x, 1, function(p)
        sprintf("(%s, %s)", p[4], p[5]))
    out <- data.frame(x)[, 1:3]
    out[["Probability interval"]] <- conf_qts
    cat(sprintf("Conjugate distribution: %s", model), sep = "\n")
    cat("", "\n")
    print(out)
    invisible()
}

#-------------------------------------------
# Probability intervals
#' @rdname methods-conjugate
#' @export
#' @param level The confidence level required.
#' @param type The type of interval.
#' @param parm Unused.
#' @importFrom rootSolve uniroot.all
# @inheritParams stats::confint
confint.conjugate <- function(object,
                              parm = "theta",
                              level = 0.95,
                              type = c("quantile", "hpd"),
                              dist = c("posterior", "likelihood",
                                       "prior"), ...) {
    type <- match.arg(type)
    dist <- match.arg(dist)
    probs <- (1 + c(-1, 1) * level) * 100 / 2
    fun <- switch(type, "quantile" = confint_qts, "hpd" = confint_hpd)
    quantis <- rbind(fun(object, level, dist))
    colnames(quantis) <- paste0(round(probs, getOption("digits")), "%")
    rownames(quantis) <- "theta"
    quantis
}
confint_qts <- function(object, level = 0.95,
                        dist = c("posterior", "likelihood",
                                 "prior")) {
    probs <- (1 + c(-1, 1) * level) / 2
    dist <- match.arg(dist)
    quantis <- object$distribution$q(probs, object[[dist]][[1]],
                                     object[[dist]][[2]])
    return(quantis)
}
confint_hpd <- function(object, level = 0.95,
                        dist = c("posterior", "likelihood",
                                 "prior")) {
    # interval to be searched for the root (I have to improve this)
    interval <- extendrange(confint_qts(object, level, dist),
                            f = 0.3)
    maxdensi <- object$distribution$d(object$summary[[dist]]$Mean,
                                      object[[dist]][[1]],
                                      object[[dist]][[2]])
    findqt <- function(x, k) {
        with(object[[dist]],
             object$distribution$d(x, alpha, beta) - k)
    }
    kchoose <- with(object[[dist]], {
        hpdfun <- function(k) {
            roots <- rootSolve::uniroot.all(findqt, interval, k = k)
            (diff(pbeta(roots, alpha, beta)) - level)^2
        }
        optimize(hpdfun, c(0, maxdensi))$minimum
    })
    quantis <- rootSolve::uniroot.all(findqt, interval, k = kchoose)
    return(quantis)
}

#-------------------------------------------
# Plots
#' @title Draw Distributions of a Univariate Bayesian Model
#' @description Draw distributions (prior, likelihood, posterior and
#'     predictive) of a univariate bayesian model.
#' @param x An object of class \code{conjugate}.
#' @param dist The distribution that is used, \code{posterior},
#'     \code{prior}, \code{likelihood} or \code{predictive} distribution.
#' @param interval interval to draw density. A list with arguments to
#'     the function \code{seq}.
#' @param type,size,ylab,xlab,... Graphics parameters.
#' @importFrom lattice xyplot
#' @export
plot.conjugate <- function(x,
                           dist = c("posterior", "likelihood",
                                    "prior", "predictive"),
                           interval = NULL,
                           type = "l",
                           size = NULL,
                           ylab = NULL,
                           xlab = NULL,
                           ...) {
    #-------------------------------------------
    # Protections
    if (!is.null(interval) & !is.list(interval)) {
        msg <- paste("`interval` should be a list with arguments",
                     "to the function `seq`")
        stop(msg)
    }
    #-------------------------------------------
    dist <- match.arg(dist)
    if (is.null(ylab)) ylab <- "Density"
    if (dist == "predictive") {
        if (is.null(xlab)) xlab <- expression(y)
        fun <- x[[dist]]
        y <- do.call(seq, interval)
        if (x$model == "Beta-Binomial") {
            fy <- fun(y, size)
        } else {
            fy <- fun(y)
        }
        out <- lattice::xyplot(fy ~ y, type = type,
                               ylab = ylab, xlab = xlab, ...)
    } else {
        if (is.null(xlab)) xlab <- expression(theta)
        pars <- x[[dist]]
        if (is.null(interval)) {
            range <- confint_qts(x, level = 1 - 1e-5, dist = dist)
            interval <- list(range[1], range[2], length = 101)
        }
        theta <- do.call(seq, interval)
        densi <- x$distribution$d(theta, pars[[1]], pars[[2]])
        out <- lattice::xyplot(densi ~ theta, type = type,
                               ylab = ylab, xlab = xlab, ...)
    }
    attr(out, "class") <- c("plot.conjugate", "trellis")
    return(out)
}

#' @title Produce a Legend or Key
#' @description Produce a legend or key based on a list of arguments.
#' @param lines,... Parameter into a list to create key.
#' @export
#' @import lattice
#' @import latticeExtra
draw_legend <- function(lines = TRUE, ...) {
    out <- list(lines = lines, ...)
    attr(out, "class") <- "plot.conjugate"
    attr(out, "islegend") <- TRUE
    return(out)
}

`+.plot.conjugate` <- function(object1, object2) {
    if (is.null(attr(object2, "islegend"))) {
        out <- latticeExtra:::`+.trellis`(object1, object2)
    } else {
        out <- lattice:::update.trellis(object1, key = object2)
    }
    return(out)
}
