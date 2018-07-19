############################################################
# OPEN INTERVAL SEARCH

search.open <- function(W0, lo=1e-4, hi=1e0, log=TRUE, ...) {
    fiter <- 2

    trycrit <- function(x) {
        y <- opencrit(x, W0=W0, log=TRUE, ...)
        if(is.na(y)) {
            y <- opencrit(x, W0=W0, log=FALSE, ...)
        }
        y
    }

    # minimum must be positive, maximum must be negative
    locont <- TRUE
    while(locont) {
        fiter <- fiter+1
        (flo <- opencrit(lo, W0=W0, log=log, ...))
        if(is.na(flo)) {
            lo <- (4/3)*lo
        } else {
            if(flo > 0) { locont <- FALSE } else { lo <- 0.5*lo }
        }
    }
    hicont <- TRUE
    while(hicont) {
        fiter <- fiter+1
        (fhi <- opencrit(hi, W0=W0, log=log, ...)) > 0
        if(is.na(fhi)) {
            hi <- (2/3)*hi
        } else {
            if(fhi < 0) { hicont <- FALSE } else { hi <- 2*hi }
        }
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
bridge.Mp.crit <- function(delta, xl, xr, Wl, Wr, Mpadj=function(...) 0, ...) {
    dbar <- xr - xl - delta
    #(Wr-Wl)/(xr-xl) + dct(delta) - dct(dbar)
    (Wr-Wl)/(xr-xl) + Mpadj(delta, dbar, ...)
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
    # which half are we in? +1 = upper, -1 = lower
    half.sign <- sign(Mp0 - mid)

    # find zero @ higher end
    f.n.vals <- f.n(upperzero.int, ...)
    if(sign(prod(f.n.vals)) < 0) {
        # if upper zero is between Mp0 (w/ numerical error!) and endpoint
        upperzero <- uniroot(f.n, upperzero.int, f.lower=f.n.vals[1], f.upper=f.n.vals[2], ...)
        res <- list(u0=upperzero$root)
    } else {
        # else, check between midpoint and Mp0
        # but! due to numerical error, f.n(Mp0) and f.f(Mp0) may not share sign
        Mp0.br <- Mp0 + Mp0-mid
        f.f.mid <- f.f(mid, ...)
        f.f.Mp0 <- f.f(Mp0, ...)
        f.f.Mp0.br <- f.f(Mp0.br, ...)

        if(half.sign > 0) {
            ends <- c(mid, Mp0)
            ends.br <- c(mid, Mp0.br)
            f.f.ends <- c(f.f.mid, f.f.Mp0)
            f.f.ends.br <- c(f.f.mid, f.f.Mp0.br)
        } else {
            ends <- c(Mp0, mid)
            ends.br <- c(Mp0.br, mid)
            f.f.ends <- c(f.f.Mp0, f.f.mid)
            f.f.ends.br <- c(f.f.Mp0.br, f.f.mid)
        }

        if(!any(is.na(f.f.ends)) && (sign(f.f.ends[1]) * sign(f.f.ends[2])) <= 0) {
            upperzero <- uniroot(f.f, ends, f.lower=f.f.ends[1], f.upper=f.f.ends[2], ...)
        #} else if (!any(is.na(f.f.ends.br)) && (sign(f.f.ends.br[1]) * sign(f.f.ends.br[2])) <= 0) {
            #upperzero <- uniroot(f.f, ends.br, f.lower=f.f.ends.br[1], f.upper=f.f.ends.br[2], ...)
        } else {
            stop('bad news')
        }
        #upperzero <- uniroot(f.f, ends, ...)
        res <- list(u0=upperzero$root)
    }

    # find zero @ lower end?
    lowermin <- optimize(f.n, loweroptim.int, ...)

    # lower min exists if objective negative in left half or positive in right half
    #if(half.sign * lowermin$objective < 0) {
    if(sign(mid - lowermin$minimum) * lowermin$objective < 0) {
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
search.brid.tight <- function(f.n, f.f, lo, mid, hi, eps=1e-6, ...) {
    midm <- f.f(mid-eps, ...)
    midp <- f.f(mid+eps, ...)
    cat(sprintf('%s %s %s %s\n', mid-eps, midm, mid+eps, midp))
    if(midm > midp) {
        upperzero <- uniroot(f.f, c(mid-eps, mid+eps),
                             f.lower=midm, f.upper=midp, ...)
        #res <- list(u0=mid)
        res <- list(u0=upperzero$root)
    } else {
        lowerzero <- uniroot(f.f, c(lo, mid-eps), ...)
        upperzero <- uniroot(f.f, c(mid+eps, hi), ...)
        res <- list(u0=upperzero$root, l0=lowerzero$root, ...)
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
    Mplo <- bridge.Mp.crit(lo, xl, xr, Wl, Wr, ...)
    Mphi <- bridge.Mp.crit(hi, xl, xr, Wl, Wr, ...)
    if(Mplo * Mphi >= 0) {
        stop('Mp has equal sign near both endpoints: probably eps too big')
    } else {
        urMp <- uniroot(bridge.Mp.crit, c(lo, hi),
                        f.lower=Mplo, f.upper=Mphi,
                        xl=xl, xr=xr, Wl=Wl, Wr=Wr, ...)
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
    #if(abs(Wr-Wl) < 1e-7){
    #if(abs(Mp0 - mid) < 1e-7){
    if(FALSE){
        res <- search.brid.tight(bridcrit.f, bridcrit.f.normal, lo, mid, hi, eps=eps, ...)
    } else {
        # which half are we in? +1 = upper, -1 = lower
        half.sign <- sign(Mp0 - mid)
        #if(Wl < Wr) {
        if(half.sign>0) {
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

        res <- tryCatch({
            res <- search.brid.regular(bridcrit.f, bridcrit.f.normal,
                                       upperzero.int, loweroptim.int, lowerzero.end,
                                       mid, Mp0, ...)
            res
        }, error=function(e) {
            # optimistically: we have guessed the wrong half, so let's try in the other?
            if(half.sign<0) {
                loweroptim.int <- c(lo, mid-eps)
                lowerzero.end <- lo
                upperzero.int <- c(Mp0, hi)
            } else {
                loweroptim.int <- c(mid+eps, hi)
                lowerzero.end <- hi
                upperzero.int <- c(lo, Mp0)
            }
            res <- search.brid.regular(bridcrit.f, bridcrit.f.normal,
                                       upperzero.int, loweroptim.int, lowerzero.end,
                                       mid, Mp0, ...)
            res
        })
        res
    }

    res
}
