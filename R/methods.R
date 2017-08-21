#' @name methods-conjugate
#' @title S3 Methods for the Object of Conjugate Class
#' @description Methods functions
#' @param x An object of class \code{conjugate}.
#' @param object An object of class \code{conjugate}.
#' @param ... unused.
#'
NULL

#' @rdname methods-conjugate
#' @author Eduardo Junior <edujrrib@gmail.com>
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
    cat(sprintf("Conjugate model: %s\n", x$model), sep = "\n")
    cat(paste0("  Prior:      ", dist, pars[[1]]), sep = "\n")
    cat(paste0("  Likelihood: ", dist, pars[[2]]), sep = "\n")
    cat(paste0("  Posterior:  ", dist, pars[[3]]), sep = "\n")
    cat("", sep = "\n")
    invisible()
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

#' @rdname methods-conjugate
#' @author Eduardo Junior <edujrrib@gmail.com>
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
    colnames(quantis) <- paste0(probs, "%")
    rownames(quantis) <- "theta"
    quantis
}
