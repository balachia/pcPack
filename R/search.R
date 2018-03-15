############################################################
# OPEN INTERVAL SEARCH

search.open <- function(W0, lo=1e-4, hi=1e5, log=TRUE, ...) {
    fiter <- 2

    # minimum must be positive, maximum must be negative
    while((flo <- opencrit(lo, W0=W0, log=log, ...)) < 0) {
        fiter <- fiter+1
        lo <- 0.5*lo
    }
    while((fhi <- opencrit(hi, W0=W0, log=log, ...)) > 0) {
        fiter <- fiter+1
        hi <- 2*hi
    }

    # search
    tryCatch({
        urres <- uniroot(opencrit, c(lo, hi),
                         f.lower=flo, f.upper=fhi,
                         W0=W0, log=log, ...)
    }, error=function(e) {
        stop(e)
    })
    res <- list(root=urres$root, preiter=fiter, uriter=urres$iter, uniroot=urres)
    res
}

############################################################
# BRIDGE INTERVAL SEARCH

#' Bridge function Mp criterion
#'
#' critical zero occurs in the Mp part of the bridge criterion
#' this helps us find it
#' @param delta distance from left endpoint
#' @param xl left endpoint
#' @param xr right endpoint
#' @param Wl value at left endpoint
#' @param Wr value at right endpoint
bridge.Mp.crit <- function(delta, xl, xr, Wl, Wr) {
    dbar <- xr - xl - delta
    (Wr-Wl)/(xr-xl) + dct(delta) - dct(dbar)
}

#' Regular search procedure for bridge interval.
#'
#' We look for maxima of the underlying function by looking for zeros of its derivative.
#' In general, there is either one maximum, or two maxima and one minimum
#' In latter case, one half of space contains the single higher valued maximum,
#' while other half of space contains a minimum near the endpoint and a maximum further away.
#' So we find maximum in upper half, and in the lower half, we look for minimum of criterion function to check if the minima/maxima crossing points exist; if they do, we find them.
#' @param f.n normal criterion function
#' @param f.f fallback criterion function
#' @param upperzero.int interval containing higher value zero
#' @param loweroptim.int interval containing function minimum in lower half
#' @param lowerzero.end endpoint of interval containing lower zero
#' @param mid midpoint
#' @param Mp0 crossing point of Mp (derivative of mean term) function component
#' @param eps unnecessary argument?
#' @param ... additional arguments to criterio function
#' @importFrom stats uniroot optimize
search.brid.regular <- function(f.n, f.f,
                                upperzero.int, loweroptim.int, lowerzero.end,
                                mid, Mp0,
                                eps=1e-6, ...) {
    # find zero @ higher end
    if(sign(prod(f.n(upperzero.int))) < 0) {
        upperzero <- uniroot(f.n, upperzero.int, ...)
        res <- list(u0=upperzero$root)
    } else {
        upperzero <- uniroot(f.f, sort(c(mid, Mp0)), ...)
        res <- list(u0=upperzero$root)
    }

    # find zero @ lower end?
    lowermin <- optimize(f.n, loweroptim.int, ...)

    # lower min exists if objective negative in left half or positive in right half
    if(sign(upperzero.int[1] - lowermin$minimum) * lowermin$objective < 0) {
        lowerzero <- uniroot(f.n, sort(c(lowerzero.end, lowermin$minimum)), ...)
        res <- c(res, list(l0=lowerzero$root))
    }

    res
}

#' Search procedure when Mp0 (mean derivative crossover) is very close to the midpoint
#'
#' bridcrit has zeros at midpoint and at mp0, so this causes problems
#' @param f criterion function
#' @param lo left endpoint of search interval
#' @param mid midpoint of search interval
#' @param hi right endpoint of search interval
#' @param eps minimum distance from endpoints to search at
#' @param ... additional arguments to criterion function (?)
#' @importFrom stats uniroot
search.brid.tight <- function(f, lo, mid, hi, eps=1e-6, ...) {
    midm <- f(mid-eps)
    midp <- f(mid+eps)
    if(midm > midp) {
        upperzero <- uniroot(f, c(mid-eps, mid+eps),
                             f.lower=midm, f.upper=midp)
        #res <- list(u0=mid)
        res <- list(u0=upperzero$root)
    } else {
        lowerzero <- uniroot(f, c(lo, mid-eps))
        upperzero <- uniroot(f, c(mid+eps, hi))
        res <- list(u0=upperzero$root, l0=lowerzero$root)
    }
    res
}

#' Search bridge interval for optima
#'
#' @param xl left endpoint
#' @param xr right endpoint
#' @param Wl value at left endpoint
#' @param Wr value at right endpoint
#' @param eps minimum distance from endpoints to search at
#' @param debug report debugging information?
#' @param ... additional arguments to criterion function
#' @importFrom stats uniroot
search.brid <- function(xl, xr, Wl, Wr, eps=1e-6, debug=FALSE, ...) {
    # for unequal Ws, we have two (three) search regions
    # on half with bigger W
    #   Mp crosses 0 somewhere
    #   we search from Mp crossover to end of region
    # on half with smaller W
    #   criterion either cross zero twice, or not at all
    #       if twice, zero closest to midpoint is minimum (no good)
    #   we search for minimum of criterion function
    #   if it is below zero, we have two crossovers,
    #       so we search from endpoint to min, and min to midpoint

    # for equal Ws, we have two search regions
    #   we either have maximum at midpoint
    #   or one maximum in each half

    # key points
    lo <- eps
    hi <- xr-xl-eps
    mid <- (xr-xl)/2

    # look for Mp crossover
    Mplo <- bridge.Mp.crit(lo, xl, xr, Wl, Wr)
    Mphi <- bridge.Mp.crit(hi, xl, xr, Wl, Wr)
    if(Mplo * Mphi >= 0) {
        stop('Mp has equal sign near both endpoints: probably eps too big')
    } else {
        urMp <- uniroot(bridge.Mp.crit, c(lo, hi),
                        f.lower=Mplo, f.upper=Mphi,
                        xl=xl, xr=xr, Wl=Wl, Wr=Wr)
        Mp0 <- urMp$root
    }

    # set up partial criteria functions
    bridcrit.f <- function(delta, ...) {
        bridcrit(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, log=TRUE, ...)
    }

    bridcrit.f.normal <- function(delta, ...) {
        bridcrit(delta, xl=xl, xr=xr, Wl=Wl, Wr=Wr, log=FALSE, ...)
    }

    # we seem to run into numerical issues once Wr-Wl gap is too small
    # empirically, about 1e-7
    if(abs(Wr-Wl) < 1e-7){
        res <- search.brid.tight(bridcrit.f, lo, mid, hi, eps=eps)
    } else {
        if(Wl < Wr) {
            # optimum in (weakly) upper half, Mp0 >= mid
            #stopifnot(Mp0 >= mid)
            loweroptim.int <- c(lo, mid-eps)
            lowerzero.end <- lo
            upperzero.int <- c(Mp0, hi)
        } else {
            # optimum in (strictly) lower half, Mp0 < mid
            #stopifnot(Mp0 < mid)
            loweroptim.int <- c(mid+eps, hi)
            lowerzero.end <- hi
            upperzero.int <- c(lo, Mp0)
        }

        if(debug) {
            cat(sprintf('%s --- %s --- %s\n', lo, mid, hi))
            cat(sprintf('Mp0: %s --- u0: %s --- lo: %s\n',
                        Mp0,
                        paste(upperzero.int, collapse=', '),
                        paste(loweroptim.int, collapse=', ')))
        }

        res <- search.brid.regular(bridcrit.f, bridcrit.f.normal,
                                   upperzero.int, loweroptim.int, lowerzero.end,
                                   mid, Mp0, ...)
    }

    res
}

# @import data.table
test_bridge <- function(xlseq=0, xrseq=10^seq(-5,5,length.out=1+1e2),
                        Wlseq=0, Wrseq=10^seq(-10,5,length.out=1+1e2),
                        ...) {
    ncol <- 50
    ngap <- 100
    params <- CJ(xl=xlseq, xr=xrseq, Wl=Wlseq, Wr=Wrseq)
    res <- sapply(1:nrow(params), function(i) {
        if(i %% ngap == 0) cat('.')
        if(i %% (ngap*ncol) == 0) cat(sprintf(' %3.2f%%\n', 100*i/nrow(params)))
        tryCatch({
            xl=params[i, xl]
            xr=params[i, xr]
            Wl=params[i, Wl]
            Wr=params[i, Wr]
            res <- search.brid(xl, xr, Wl, Wr, ...)
            ''
        }, error=function(e) {
            #cat(sprintf('%s, %s: %s, %s: ', xl, xr, Wl, Wr))
            #cat('\n')
            #print(e$message)
            #print(e$call)
            e$message
        })
    })
    cat('\n')
    params[, res := res]
    params
}

bridge_test_diagnostic <- function(dt) {
    ggp <- ggplot(dt, aes(xr, Wr, fill=res)) +
        theme(legend.position='bottom') +
        scale_x_log10() +
        scale_y_log10() +
        coord_cartesian(expand=FALSE) +
        scale_fill_brewer(palette='Dark2') +
        geom_tile(alpha=0.75)
    ggp
}
