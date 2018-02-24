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

#' Expected first derivative of utility function
#'
#' risk aversion = b * (1 / (1 + a/b exp(bm)))
#' a controls straightness, b controls bendiness
#' RA increasing in b, decreasing in a
#' @param m wealth
#' @param a a parameter
#' @param b b parameter
Edu <- function(delta, W0, a=1, b=1, ...) {
    a + b * exp(-b * (W0 + ct(delta)) + 0.5 * b^2 * sigma^2 * delta)
}

#' Expected second derivative of utility function
#'
#' risk aversion = b * (1 / (1 + a/b exp(bm)))
#' a controls straightness, b controls bendiness
#' RA increasing in b, decreasing in a
#' @param m wealth
#' @param a a parameter
#' @param b b parameter
Ed2u <- function(delta, W0, a=1, b=1, ...) -b^2 * exp(-b * (W0 + ct(delta)) + 0.5 * b^2 * sigma^2 * delta)

#' Expected utility on open interval
#'
#' @param delta distance from known point
#' @param W0 value at known point
#' @param a utility function *a* parameter
#' @param b utility function *b* parameter
Euopen <- function(delta, W0, sigma=1, a=1, b=1, ...) {
    exparg <- exparg.open(delta, W0=W0, sigma=sigma, b=b)
    #a * (W0 + ct(delta)) - exp(-b * (W0 + c(delta)) + 0.5 * b^2 * sigma^2 * delta)
    a * (W0 + ct(delta)) - exp(exparg)
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
Eubrid <- function(delta, xl, xr, Wl, Wr, sigma=1, a=1, b=1, ...) {
    exparg <- exparg.brid(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b)
    dbar <- xr - xl - delta
    arg <- delta * (Wr - Wl) / (xr - xl) + Wl + ct(delta) + ct(dbar)
    a * arg - exp(exparg)
}

#' New value at open jump
#'
#' @param delta jump distance
#' @param sigma brownian standard deviation
#' @export
openjump <- function(delta, sigma=1) sigma * sqrt(delta) * rnorm(1)

#' New value at bridge jump
#'
#' @param delta jump distance from left endpoint
#' @param xl left endpoint
#' @param xl right endpoint
#' @param Wl value of left point
#' @param Wr value of right point
#' @param sigma brownian standard deviation
#' @export
bridjump <- function(delta, xl, xr, Wl, Wr, sigma=1) Wl + delta * (Wr - Wl) / (Wr - Wl) + sigma * sqrt((delta * (xr - xl - delta)) / (xr - xl)) * rnorm(1)



#' Ratio of -Ed2u / Edu
#' 
#' Multiplier on variance term, used in criteria calculations.
#' Logarithmic version relies on log(1 + exp(x)) for positive x.
#' For large x, this ~= x + exp(-x)
#' @param exparg argument inside exponential term
#' @param a utility function *a* parameter
#' @param b utility function *b* parameter
#' @param log use logarithmic version?
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
euratio.normal <- function(exparg, a=1, b=1, ...) {
    b^2 / a * exp(exparg) / (1 + (b/a) * exp(exparg))
}

# in general, for mean M and variance V,
# exparg = -b * M + (1/2) * b^2 * V

# in general, criterion is 0 = M' + ( 1/2 * V' * ( -Ed2u/Edu ) )
# logged criterion is 0 = -log M' + log(1/2 V') + log( -Ed2u/Edu )
# we denote these as nlogMp + loghVp + logEurat()

exparg.open <- function(delta, W0, sigma=1, b=1) {
    exparg <- -b * (W0 + ct(delta)) + 0.5 * b^2 * sigma^2 * delta
    exparg
}

exparg.brid <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, ...) {
    dbar <- xr - xl - delta
    exparg <- -b * (Wl + ct(delta) + ct(dbar) + delta*((Wr-Wl)/(xr-xl))) + 
        0.5 * b^2 * sigma^2 * ((dbar * delta) / (xr - xl))
    exparg
}

opencrit.normal <- function(delta, W0, sigma=1, b=1, debug=FALSE, ...) {
    exparg <- exparg.open(delta, W0=W0, sigma=sigma, b=b)
    Mp <- dct(delta)
    hVp <- 0.5 * sigma^2
    eurat <- euratio.normal(exparg, b=b, ...)
    crit <- Mp - hVp*eurat
    if(debug) {
        crit <- list(crit=crit, Mp=Mp, hVp=hVp, eurat=eurat, hVpeurat=hVp*eurat)
    }
    crit
}

opencrit.logged <- function(delta, W0, sigma=1, b=1, debug=FALSE, ...) {
    exparg <- exparg.open(delta, W0=W0, sigma=sigma, b=b)
    logMp <- dct.logged(delta)
    loghVp <- 2*log(sigma) - log(2)
    eurat <- euratio.logged(exparg, b=b, ...)
    crit <- logMp - loghVp - eurat
    if(debug) {
        crit <- list(crit=crit, Mp=logMp, hVp=loghVp, eurat=eurat, hVpeurat=loghVp+eurat)
    }
    crit
}

bridcrit.normal <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, debug=FALSE, ...) {
    dbar <- xr - xl - delta
    exparg <- exparg.brid(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    Mp <- ((Wr-Wl)/(xr-xl)) + dct(delta) - dct(dbar)
    hVp <- 0.5 * sigma^2 * (dbar-delta)/(xr-xl)
    eurat <- euratio.normal(exparg, b=b, ...)
    crit <- Mp - hVp*eurat
    if(debug) {
        crit <- list(crit=crit, Mp=Mp, hVp=hVp, eurat=eurat, hVpeurat=hVp*eurat)
    }
    crit
}

bridcrit.logged <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, debug=FALSE, ...) {
    dbar <- xr - xl - delta
    # flip signs in back half of search space
    flip <- sign(dbar-delta)
    exparg <- exparg.brid(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    Mp <- ((Wr-Wl)/(xr-xl)) + dct(delta) - dct(dbar)
    logMp <- ifelse(flip*Mp>0, log(flip*Mp), -Inf)
    loghVp <- 2*log(sigma) - log(2) - log(xr-xl) + log(flip*(dbar-delta))
    eurat <- euratio.logged(exparg, b=b, ...)
    logMp <- flip*logMp
    loghVp <- flip*loghVp
    eurat <- flip*eurat
    crit <- logMp - loghVp - eurat
    #crit <- flip*(logMp - loghVp - eurat)
    if(debug) {
        crit <- list(crit=crit, Mp=logMp, hVp=loghVp, eurat=eurat, hVpeurat=loghVp+eurat)
    }
    crit
}

#' open interval criterion function
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
opencrit <- function(delta, W0, sigma=1, b=1, log=TRUE, ...) {
    if(log) {
        opencrit.logged(delta, W0=W0, sigma=sigma, b=b, ...)
    } else {
        opencrit.normal(delta, W0=W0, sigma=sigma, b=b, ...)
    }
}

#' bridge interval criterion function
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
bridcrit <- function(delta, xl, xr, Wl, Wr, sigma=1, b=1, log=TRUE, ...) {
    if(log) {
        bridcrit.logged(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    } else {
        bridcrit.normal(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, sigma=sigma, b=b, ...)
    }
}

# @import data.table
plot_crit <- function(critfs, n=1e3+1, ylim=0.95, xl=0, xr=1, ...) {
    if(!'list' %in% class(critfs)) critfs <- list(critfs)
    xs <- seq(xl, xr, length.out=n)

    ysdts <- list()
    for(i in seq_along(critfs)) {
        critf <- critfs[[i]]
        ys <- critf(xs, xl=xl, xr=xr, debug=TRUE, ...)
        cols <- c('crit', 'Mp', 'hVp', 'eurat', 'hVpeurat')
        ysdt <- rbindlist(map(cols, ~ data.table(x=xs, y=ys[[.]], type=., func=i)))
        ysdts <- c(ysdts, list(ysdt))
    }
    ysdt <- rbindlist(ysdts)

    ysdt[, maxy := quantile(abs(y), ylim, na.rm=TRUE), by=func]
    plotdt <- ysdt[abs(y) <= maxy | is.na(y) | is.infinite(y)]

    plotcols <- cols[c(2,5)]
    ggp <- ggplot(plotdt[type %in% plotcols], aes(x, y, color=type)) +
        facet_wrap(~func, ncol=1, scales='free_y') +
        geom_hline(yintercept=0) +
        #coord_cartesian(ylim=c(-ylim, ylim)) +
        coord_cartesian(expand=FALSE) +
        geom_line() +
        geom_line(data=plotdt[type=='crit'], color='black', linetype='22')
    print(ggp)

    ysdt
}

plot_bridcrit <- function(xl, xr, Wl, Wr, ...) {
    plot_crit(bridcrit, xl=xl, xr=xr, Wl=Wl, Wr=Wr, ...)
}

plot_bridcrits <- function(xl, xr, Wl, Wr, ...) {
    plot_crit(list(bridcrit.normal, bridcrit.logged),
              xl=xl, xr=xr, Wl=Wl, Wr=Wr, ...)
}
