#' @author Eduardo Junior <edujrrib@gmail.com>
#' @title Parser of the Model Name to Probability Distribution
#' @description Get the probability distribution functions based on
#'     conjugacy name.
#' @param model An object of class \code{conjugate}.
#' @return The [d, p, q, r] functions of the congate distribution.
#' @importFrom grDevices extendrange
#' @importFrom stats dbeta dgamma dnorm pbeta pgamma pnorm qbeta qgamma
#'     qnorm rbeta rgamma rnorm
#'
model2dist <- function(model) {
    out <- switch(
        model,
        "Beta-Binomial"  = list("d" = dbeta,
                                "p" = pbeta,
                                "q" = qbeta,
                                "r" = rbeta),
        "Gamma-Poisson"  = list("d" = dgamma,
                                "p" = pgamma,
                                "q" = qgamma,
                                "r" = rgamma),
        "Normal-Normal"  = list("d" = dnorm,
                                "p" = pnorm,
                                "q" = qnorm,
                                "r" = rnorm),
        "Gamma-Gamma"    = list("d" = dgamma,
                                "p" = pgamma,
                                "q" = qgamma,
                                "r" = rgamma),
        "Beta-Geometric" = list("d" = dbeta,
                                "p" = pbeta,
                                "q" = qbeta,
                                "r" = rbeta),
        "BinomNeg-Beta"  = list("d" = dbeta,
                                "p" = pbeta,
                                "q" = qbeta,
                                "r" = rbeta)
    )
    out
}
