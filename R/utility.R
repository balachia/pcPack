############################################################
# COMPETITION 

#' Competition function
#'
#' Competition function
#' @param delta distance from competitor
ct <- function(delta) -1 / delta

#' Competition function (d/ddelta)
#'
#' Competition function (d/ddelta)
#' @param delta distance from competitor
dct <- function(delta) 1 / (delta^2)

#' Log competition function (d/ddelta)
#'
#' Log competition function (d/ddelta)
#' @param delta distance from competitor
dct.logged <- function(delta) -2*log(delta)

ct.Madj.open <- function(delta, ...) ct(delta)
ct.Mpadj.open <- function(delta, ...) dct(delta)
ct.Madj.brid <- function(delta, dbar, ...) ct(delta) + ct(dbar)
ct.Mpadj.brid <- function(delta, dbar, ...) dct(delta) - dct(dbar)

############################################################
# Helpers

#' Squash all arguments to 0
#'
squash <- function(...) 0

############################################################
# UTILITY

#' Utility function
#'
#' risk aversion = b * (1 / (1 + a/b exp(bm)))
#' a controls straightness, b controls bendiness
#' RA increasing in b, decreasing in a
#' @param m wealth
#' @param a a parameter
#' @param b b parameter
u <- function(m, a=1, b=1) a * m - exp(-b * m)

# #' Expected first derivative of utility function
# #'
# #' risk aversion = b * (1 / (1 + a/b exp(bm)))
# #' a controls straightness, b controls bendiness
# #' RA increasing in b, decreasing in a
# #' @param m wealth
# #' @param a a parameter
# #' @param b b parameter
# #' @param sigma random walk sigma parameter
# Edu <- function(delta, W0, a=1, b=1, sigma=1, ...) {
#     a + b * exp(-b * (W0 + ct(delta)) + 0.5 * b^2 * sigma^2 * delta)
# }

# #' Expected second derivative of utility function
# #'
# #' risk aversion = b * (1 / (1 + a/b exp(bm)))
# #' a controls straightness, b controls bendiness
# #' RA increasing in b, decreasing in a
# #' @param m wealth
# #' @param a a parameter
# #' @param b b parameter
# #' @param sigma random walk sigma parameter
# Ed2u <- function(delta, W0, a=1, b=1, ...) -b^2 * exp(-b * (W0 + ct(delta)) + 0.5 * b^2 * sigma^2 * delta)

#' Expected utility on open interval
#'
#' @param delta distance from known point
#' @param W0 value at known point
#' @param sigma random walk sigma parameter
#' @param a utility function *a* parameter
#' @param b utility function *b* parameter
#' @param ... capture additional arguments
Euopen <- function(delta, W0, sigma=1, a=1, b=1, Madj=squash, ...) {
    #exparg <- exparg.open(delta, W0=W0, sigma=sigma, b=b)
    exparg <- exparg.open(delta, W0=W0, sigma=sigma, b=b, Madj=Madj, ...)
    #a * (W0 + ct(delta)) - exp(-b * (W0 + c(delta)) + 0.5 * b^2 * sigma^2 * delta)
    #a * (W0 + ct(delta)) - exp(exparg)
    a * (W0 + Madj(delta, ...)) - exp(exparg)
}

#' Expected utility on bridge interval
#'
#' @param delta distance from left point
#' @param xl location of left point
#' @param xr location of right point
#' @param Wl value of left point
#' @param Wr value of right point
#' @param sigma brownian walk *sigma* parameter
#' @param a utility function *a* parameter
#' @param b utility function *b* parameter
#' @param ... capture additional arguments
Eubrid <- function(delta, xl, xr, Wl, Wr, sigma=1, a=1, b=1, Madj=squash, ...) {
    #exparg <- exparg.brid(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b)
    exparg <- exparg.brid(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, Madj=Madj, ...)
    dbar <- xr - xl - delta
    #arg <- delta * (Wr - Wl) / (xr - xl) + Wl + ct(delta) + ct(dbar)
    arg <- delta * (Wr - Wl) / (xr - xl) + Wl + Madj(delta, dbar, ...)
    a * arg - exp(exparg)
}

#' New value at open jump
#'
#' @param delta jump distance
#' @param sigma brownian standard deviation
#' @export
#' @importFrom stats rnorm
openjump <- function(delta, sigma=1, ...) sigma * sqrt(delta) * rnorm(1)

#' New value at bridge jump
#'
#' @param delta jump distance from left endpoint
#' @param xl left endpoint
#' @param xr right endpoint
#' @param Wl value of left point
#' @param Wr value of right point
#' @param sigma brownian standard deviation
#' @export
#' @importFrom stats rnorm
bridjump <- function(delta, xl, xr, Wl, Wr, sigma=1, ...) Wl + delta * (Wr - Wl) / (xr - xl) + sigma * sqrt((delta * (xr - xl - delta)) / (xr - xl)) * rnorm(1)



#' Ratio of -Ed2u / Edu
#' 
#' Multiplier on variance term, used in criteria calculations.
#' Logarithmic version relies on log(1 + exp(x)) for positive x.
#' For large x, this ~= x + exp(-x)
#' @param exparg argument inside exponential term
#' @param a utility function *a* parameter
#' @param b utility function *b* parameter
#' @param log use logarithmic version?
#' @param ... additional arguments to subfunctions
#' @export
euratio <- function(exparg, a=1, b=1, log=TRUE, ...) {
    if(log) {
        euratio.logged(exparg, a=a, b=b, ...)
    } else {
        euratio.normal(exparg, a=a, b=b, ...)
    }
}

#' Logarithm of -Ed2u / Edu
#' 
#' Relies on log(1 + exp(x)) for positive x. For large x, this ~= x + exp(-x)
#' @param exparg argument inside exponential term
#' @param a utility function *a* parameter
#' @param b utility function *b* parameter
#' @param ... capture additional arguments
euratio.logged <- function(exparg, a=1, b=1, ...) {
    exparg2 <- exparg + log(b/a)
    #log1pexp <- if(exparg2 <= 18) { log1p(exp(exparg2)) } else { exparg2 + exp(-exparg2) }
    #log1pexp <- log1p(exp(exparg2))*(exparg <= 18) + (exparg2+exp(-exparg2))*(exparg>18)
    log1pexp <- ifelse(exparg2 <= 18, log1p(exp(exparg2)), exparg2 + exp(-exparg2))

    2*log(b) - log(a) + exparg - log1pexp
}

#' -Ed2u / Edu
#' 
#' @param exparg argument inside exponential term
#' @param a utility function *a* parameter
#' @param b utility function *b* parameter
#' @param ... capture additional arguments
euratio.normal <- function(exparg, a=1, b=1, ...) {
    b^2 / a * exp(exparg) / (1 + (b/a) * exp(exparg))
}

# in general, for mean M and variance V,
# exparg = -b * M + (1/2) * b^2 * V

# in general, criterion is 0 = M' + ( 1/2 * V' * ( -Ed2u/Edu ) )
# logged criterion is 0 = -log M' + log(1/2 V') + log( -Ed2u/Edu )
# we denote these as nlogMp + loghVp + logEurat()

#' Calculate Eu exponential argument on open interval
#'
exparg.open <- function(delta, W0, sigma=1, b=1, Madj=squash, ...) {
    #exparg <- -b * (W0 + ct(delta)) + 0.5 * b^2 * sigma^2 * delta
    exparg <- -b * (W0 + Madj(delta, ...)) + 0.5 * b^2 * sigma^2 * delta
    exparg
}

#' Calculate Eu exponential argument on bridge interval
#'
exparg.brid <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, Madj=squash, ...) {
    dbar <- xr - xl - delta
    #exparg <- -b * (Wl + ct(delta) + ct(dbar) + delta*((Wr-Wl)/(xr-xl))) + 
    exparg <- -b * (Wl + Madj(delta, dbar, ...) + delta*((Wr-Wl)/(xr-xl))) + 
        0.5 * b^2 * sigma^2 * ((dbar * delta) / (xr - xl))
    exparg
}

#' Calculate open interval criterion
#'
opencrit.normal <- function(delta, W0, sigma=1, b=1, Mpadj=squash, augment=FALSE, ...) {
    exparg <- exparg.open(delta, W0=W0, sigma=sigma, b=b, ...)
    #Mp <- dct(delta)
    Mp <- Mpadj(delta, ...)
    hVp <- 0.5 * sigma^2
    eurat <- euratio.normal(exparg, b=b, ...)
    crit <- Mp - hVp*eurat
    if(augment) {
        attr(crit, 'crit') <- crit
        attr(crit, 'Mp') <- Mp
        attr(crit, 'hVp') <- hVp
        attr(crit, 'eurat') <- eurat
        attr(crit, 'hVpeurat') <- hVp*eurat
    }
    crit
}

#' Calculate open interval criterion (logged)
#'
opencrit.logged <- function(delta, W0, sigma=1, b=1, Mpadj=squash, augment=FALSE, ...) {
    exparg <- exparg.open(delta, W0=W0, sigma=sigma, b=b, ...)
    #logMp <- dct.logged(delta)
    logMp <- log(Mpadj(delta, ...))
    loghVp <- 2*log(sigma) - log(2)
    eurat <- euratio.logged(exparg, b=b, ...)
    crit <- logMp - loghVp - eurat
    if(augment) {
        attr(crit, 'crit') <- crit
        attr(crit, 'Mp') <- logMp
        attr(crit, 'hVp') <- loghVp
        attr(crit, 'eurat') <- eurat
        attr(crit, 'hVpeurat') <- loghVp+eurat
    }
    crit
}

#' Calculate bridge interval criterion
#'
bridcrit.normal <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, Mpadj=squash, augment=FALSE, ...) {
    dbar <- xr - xl - delta
    exparg <- exparg.brid(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    #Mp <- ((Wr-Wl)/(xr-xl)) + dct(delta) - dct(dbar)
    Mp <- ((Wr-Wl)/(xr-xl)) + Mpadj(delta, dbar, ...)
    hVp <- 0.5 * sigma^2 * (dbar-delta)/(xr-xl)
    eurat <- euratio.normal(exparg, b=b, ...)
    crit <- Mp - hVp*eurat
    if(augment) {
        attr(crit, 'crit') <- crit
        attr(crit, 'Mp') <- Mp
        attr(crit, 'hVp') <- hVp
        attr(crit, 'eurat') <- eurat
        attr(crit, 'hVpeurat') <- hVp*eurat
    }
    crit
}

#' Calculate bridge interval criterion (logged)
#'
bridcrit.logged <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, Mpadj=squash, augment=FALSE, ...) {
    dbar <- xr - xl - delta
    # flip signs in back half of search space
    flip <- sign(dbar-delta)
    exparg <- exparg.brid(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    #Mp <- ((Wr-Wl)/(xr-xl)) + dct(delta) - dct(dbar)
    Mp <- ((Wr-Wl)/(xr-xl)) + Mpadj(delta, dbar, ...)
    #logMp <- ifelse(flip*Mp>0, log(flip*Mp), -Inf)
    # avoid warning about NaN
    fMp <- flip*Mp
    logMp <- rep(-Inf, length.out=length(Mp))
    logMp[fMp>0] <- log(fMp[fMp>0])
    loghVp <- 2*log(sigma) - log(2) - log(xr-xl) + log(flip*(dbar-delta))
    eurat <- euratio.logged(exparg, b=b, ...)
    logMp <- flip*logMp
    loghVp <- flip*loghVp
    eurat <- flip*eurat
    crit <- logMp - loghVp - eurat
    #crit <- flip*(logMp - loghVp - eurat)
    if(augment) {
        attr(crit, 'crit') <- crit
        attr(crit, 'Mp') <- logMp
        attr(crit, 'hVp') <- loghVp
        attr(crit, 'eurat') <- eurat
        attr(crit, 'hVpeurat') <- loghVp+eurat
    }
    crit
}

#' Open interval criterion function
#'
#' General criterion is 0 = M' + ( 1/2 * V' * ( -Ed2u/Edu ) )
#' Logged criterion is 0 = -log M' + log(1/2 V') + log( -Ed2u/Edu )
#' General criterion: 0 = -log M' + log(1/2 V') + log( -Ed2u/Edu )
#' Open criterion: M' = c'(delta)
#' Open criterion: V' = sigma^2
#' @param delta distance from point
#' @param W0 value at point
#' @param sigma brownian walk *sigma* parameter
#' @param b utility function *b* parameter
#' @param log use logged criterion?
#' @param ... additional arguments to criterion subfunctions
opencrit <- function(delta, W0, sigma=1, b=1, log=TRUE, ...) {
    if(log) {
        opencrit.logged(delta, W0=W0, sigma=sigma, b=b, ...)
    } else {
        opencrit.normal(delta, W0=W0, sigma=sigma, b=b, ...)
    }
}

#' Bridge interval criterion function
#'
#' General criterion is 0 = M' + ( 1/2 * V' * ( -Ed2u/Edu ) )
#' Logged criterion is 0 = -log M' + log(1/2 V') + log( -Ed2u/Edu )
#' Bridge criterion: M' = (Wr-Wl)/(xr-xl) + c'(delta) - c'(dbar)
#' Bridge criterion: V' = (dbar-delta)/(xr-xl) * sigma^2
#' @param delta distance from left point
#' @param xl position of left point
#' @param xr position of right point
#' @param Wl value at left point
#' @param Wr value at right point
#' @param sigma brownian walk *sigma* parameter
#' @param b utility function *b* parameter
#' @param log use logged criterion?
#' @param ... additional arguments to criterion subfunctions
bridcrit <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, log=TRUE, ...) {
    if(log) {
        bridcrit.logged(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    } else {
        bridcrit.normal(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    }
}
