#' @name conjugate-models
#' @author Eduardo Junior <edujrrib@gmail.com>
#' @title Bayesian Models of the Conjugate Distributions
#' @description Compute components for the bayesian models of the
#'     conjugate distributions.
#' @param prior A parameters list of the prior distribution.
#' @param data A vector of observed data or a named list with observed
#'     data (name \code{y}) and know parameters. Beta-Binomial know
#'     parameter is size (named \code{size}) and for Normal-Normal know
#'     parameter is sigma (named \code{sigma}).
#' @return An object of class \code{conjugate}.
#'
NULL

#-------------------------------------------
# Conjugacy Beta prior and Binomial likelihood (with know size)
#' @rdname conjugate-models
#' @export
#'
beta_binomial <- function(prior, data) {
    # Organize data
    if (is.list(data)) {
        y <- data[["y"]]
        m <- data[["size"]]
    } else {
        y <- data
        m <- 1
    }
    # Compute useful statistics
    n <- length(y)
    sy <- sum(y)
    if (length(m) == 1) {
        sm <- n * m
    } else {
        sm <- sum(m)
    }
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
    # Summarise distributions
    dnames <- c("prior", "likelihood", "posterior")
    names(dnames) <- dnames
    summary <- lapply(dnames, function(d) {
        with(eval(parse(text = d)), {
            mu <- alpha / (alpha + beta)
            va <- alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
            mo <- (alpha - 1) / (alpha + beta - 2)
            mo <- mo * ifelse(alpha > 1 & beta > 1, 1, NA)
            list("Mean" = mu, "Mode" = mo, "Variance" = va)
        })
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

#-------------------------------------------
# Conjugacy Gamma prior and Poisson likelihood
#' @rdname conjugate-models
#' @export
#'
gamma_poisson <- function(prior, data) {
    # Compute useful statistics
    y <- data
    sy <- sum(y)
    n <- length(y)
    # Parameters
    likelihood <- list("shape" = sy + 1, "rate" = n)
    posterior <- list("shape" = prior$shape + sy,
                      "rate" = prior$rate + n)
    # Predictive distribution
    predictive <- function(x, log = FALSE) {
        out <- with(posterior,  {
            shape * log(rate) - lgamma(shape) - lfactorial(x) +
                lgamma(x + shape) - (x + shape) * log(1 + rate)
        })
        if (!log) {
            out <- exp(out)
        }
        return(out)
    }
    # Summarise distributions
    dnames <- c("prior", "likelihood", "posterior")
    names(dnames) <- dnames
    summary <- lapply(dnames, function(d) {
        with(eval(parse(text = d)), {
            mu <- shape / rate
            va <- shape / rate^2
            mo <- (shape - 1) / rate
            mo <- mo * ifelse(shape > 1, 1, NA)
            list("Mean" = mu, "Mode" = mo, "Variance" = va)
        })
    })
    # Output
    out <- list("model" = "Gamma-Poisson",
                "distribution" = model2dist("Gamma-Poisson"),
                "likelihood" = likelihood,
                "posterior" = posterior,
                "prior" = prior,
                "predictive" = predictive,
                "summary" = summary)
    attr(out, "class") <- "conjugate"
    return(out)
}

#-------------------------------------------
# Conjugacy Normal prior and Normal likelihood (with know sigma)
#' @rdname conjugate-models
#' @export
#' @importFrom stats sd
normal_normal <- function(prior, data) {
    # Organize data
    if (is.list(data)) {
        y <- data[["y"]]
        sigma <- data[["sigma"]]
    } else {
        y <- data
        sigma <- sd(sigma)
        msg <- sprintf("Sigma parameter take to be %.3f", sigma)
    }
    # Compute useful statistics
    n <- length(y)
    sy <- sum(y)
    # Parameters
    likelihood <- list("mean" = sy / n, "sd" = sigma / sqrt(n))
    posterior <- list(
        "mean" = ((prior$mean / prior$sd^2 + sy / sigma^2) /
                  (1 / prior$sd^2 + n / sigma^2)),
        "sd" = (1 / prior$sd^2 + n / sigma^2)^(-1)
    )
    # Predictive distribution
    predictive <- function(x, log = FALSE) {
        out <- with(posterior,  {
            dnorm(x, mean, sqrt(sigma^2 + sd^2), log = log)
        })
        return(out)
    }
    # Summarise distributions
    dnames <- c("prior", "likelihood", "posterior")
    names(dnames) <- dnames
    summary <- lapply(dnames, function(d) {
        with(eval(parse(text = d)), {
            mu <- mean
            va <- sd^2
            mo <- mean
            list("Mean" = mu, "Mode" = mo, "Variance" = va)
        })
    })
    # Output
    out <- list("model" = "Normal-Normal",
                "distribution" = model2dist("Normal-Normal"),
                "likelihood" = likelihood,
                "posterior" = posterior,
                "prior" = prior,
                "predictive" = predictive,
                "summary" = summary)
    attr(out, "class") <- "conjugate"
    return(out)
}
