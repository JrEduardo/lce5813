#-------------------------------------------
# Documentation of methods
#' @name methods-conjugate
#' @title S3 Methods for the Object of Conjugate Class
#' @description Methods functions
#' @param x An object of class \code{conjugate}.
#' @param object An object of class \code{conjugate}.
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
# Summary posterior
#' @rdname methods-conjugate
#' @export
summary.conjugate <- function(object, ...) {
    probs <- c(0.025, 0.975)
    cname <- unlist(lapply(c("QTS", "HPD"), paste0, probs))
    conf_qts <- confint_qts(object, level = 0.95)
    conf_hpd <- confint_hpd(object, level = 0.95)
    out <- rbind(c(do.call(c, object$summary), conf_qts, conf_hpd))
    colnames(out)[-(1:3)] <- cname
    rownames(out) <- "theta"
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
    nc <- max(nchar(as.integer(unlist(x))))
    x <- formatC(x, digits = digits, width = nc + digits + 1,
                 format = "f")
    conf_qts <- sprintf("(%s, %s)", x[4], x[5])
    out <- data.frame(x)[, 1:3]
    out[["Confidence interval"]] <- conf_qts
    cat(sprintf("Conjugate distribution: %s", model), sep = "\n")
    cat("Posterior summary\n", sep = "\n")
    print(out)
    invisible()
}

#-------------------------------------------
# Conficence intervals
#' @rdname methods-conjugate
#' @export
#' @param level The confidence level required.
#' @param type The type of interval.
#' @param parm Unused.
#' @importFrom rootSolve uniroot.all
# @inheritParams stats::confint
confint.conjugate <- function(object, parm = "theta", level = 0.95,
                              type = c("quantile", "hpd"), ...) {
    probs <- (1 + c(-1, 1) * level) * 100 / 2
    type <- match.arg(type)
    fun <- switch(type, "quantile" = confint_qts, "hpd" = confint_hpd)
    quantis <- rbind(fun(object, level))
    colnames(quantis) <- paste0(round(probs, getOption("digits")), "%")
    rownames(quantis) <- "theta"
    quantis
}
confint_qts <- function(object, level = 0.95) {
    probs <- (1 + c(-1, 1) * level) / 2
    quantis <- with(object, {
        distribution$q(probs, posterior[[1]], posterior[[2]])
    })
    return(quantis)
}
confint_hpd <- function(object, level = 0.95) {
    # interval to be searched for the root (I have to improve this)
    interval <- extendrange(confint_qts(object, level), f = 0.3)
    maxdensi <- with(object, {
        distribution$d(summary$Mean, posterior[[1]], posterior[[2]])
    })
    findqt <- function(x, k) {
        with(object$posterior,
             object$distribution$d(x, alpha, beta) - k)
    }
    kchoose <- with(object$posterior, {
        hpdfun <- function(k) {
            roots <- rootSolve::uniroot.all(findqt, interval, k = k)
            (diff(pbeta(roots, alpha, beta)) - level)^2
        }
        optimize(hpdfun, c(0, maxdensi))$minimum
    })
    quantis <- rootSolve::uniroot.all(findqt, interval, k = kchoose)
    return(quantis)
}
