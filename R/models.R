#' @name conjugate-models
#' @author Eduardo Junior <edujrrib@gmail.com>
#' @title Bayesian Models of the Conjugate Distributions
#' @description Compute components for the bayesian models of the
#'     conjugate distributions.
#' @param prior A parameters list of the prior distribution.
#' @param y A vector (or matrix) of observed data.
#' @return An object of class \code{conjugate}.
#'
NULL

#' @rdname conjugate-models
#' @author Eduardo Junior <edujrrib@gmail.com>
#' @title Conjucacy Beta-Binomial
#' @export
#'
beta_binomial <- function(prior, y) {
    # Organize data
    if (is.null(dim(y))) {
        m <- rep(1, length(y))
    } else {
        m <- y[, 2]
        y <- y[, 1]
    }
    # Compute useful statistics
    sy <- sum(y)
    sm <- sum(m)
    n <- length(y)
    # Parameters
    likelihood <- list("alpha" = sy + 1, "beta" = sm - sy + 1)
    posterior <- list("alpha" = prior$alpha + sy,
                      "beta" = prior$beta + sm - sy)
    # Predictive distribution
    predictive <- function(x, size, log = FALSE) {
        out <- with(posterior,  {
            lchoose(size, x) +
                lbeta(alpha + x, size + beta - x) -
                lbeta(alpha, beta)
        })
        if (!log) {
            out <- exp(out)
        }
        return(out)
    }
    # Summarise posterior
    summary <- with(posterior, {
        mu <- alpha / (alpha + beta)
        va <- alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
        mo <- (alpha - 1) / (alpha + beta - 2)
        mo <- mo * ifelse(alpha > 1 & beta > 1, 1, NA)
        list("Mean" = mu, "Mode" = mo, "Variance" = va)
    })
    # Output
    out <- list("model" = "Beta-Binomial",
                "distribution" = model2dist("Beta-Binomial"),
                "likelihood" = likelihood,
                "posterior" = posterior,
                "prior" = prior,
                "predictive" = predictive,
                "summary" = summary)
    attr(out, "class") <- "conjugate"
    return(out)
}
