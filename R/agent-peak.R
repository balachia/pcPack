############################################################
# single peak deterministic agent

#' Determiknistic agent
#'
#' Deterministic agent that descends down predetermined slopes from known peaks
#' @param n number of entries into the market
#' @param slope slope from known peaks (should be negative)
#' @param ... capture additional arguments
#' @export
make_peak_agent <- function(peaks = 0, slope = -1, ...) {
    ag <- list()
    class(ag) <- c('peak.agent', 'agent')

    # slope must be less than 0 else this fails
    stopifnot(slope < 0)
    ag$parameters <- list(peaks = peaks, slope = slope)
    ag$slope <- slope
    ag$peaks <- peaks
    #ag$plan <- create_plans_table(n+1)

    ag$Madj.open <- ct.Madj.open
    ag$Mpadj.open <- ct.Mpadj.open
    ag$Madj.brid <- ct.Madj.brid
    ag$Mpadj.brid <- ct.Mpadj.brid

    ag
}

#' @export
Ops.peak.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'peak.agent') && inherits(a2, 'peak.agent')) {
                   all(a1$peaks == a2$peaks) && a1$slope == a2$slope
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.peak.agent <- function(ag, ...) {
    cat(sprintf('peak agent :: peaks %s :: slope %s\n', paste(ag$peaks, sep = ','), ag$slope))
}

agentSetup.peak.agent <- function(agent, n, ...) {
    agent$plan <- create_plans_table(n+1)
    agent
}

#' @import data.table
agentEntry.peak.agent <- function(agent, i, intervals, positions, ...) {
    enteridx <- agent$plan[which.max(Eu), id]
    delta <- agent$plan[enteridx, delta]

    # entering at open or closed interval?
    xl <- intervals[enteridx, xl]
    xr <- intervals[enteridx, xr]
    if(is.infinite(xl) || is.infinite(xr)) {
        if(is.infinite(xl)) {
            x <- xr-delta
            W <- agent$slope * (agent$peak - x)
        } else {
            x <- xl+delta
            W <- agent$slope * (x - agent$peak)
        }
    } else {
        x <- xl+delta
        W <- agent$slope * abs(x - agent$peak)
    }

    list(idx=enteridx, x=x, W=W)
}

#' @import data.table
agentUpdate.peak.agent <- function(agent, update.idx, intervals, ...) {
    search.open.peak <- function() sqrt(-1 / agent$slope)
    search.brid.peak.f <- function(delta, xl, xr, dir = 1) dir * agent$slope * delta + ct(delta) + ct(xr - xl - delta)
    search.brid.peak <- function(xl, xr, dir = 1) optimize(search.brid.peak.f, xl = xl, xr = xr, dir = dir,
                                                             lower = 0, upper = xr - xl, maximum = TRUE)
    plan <- agent$plan
    for(idx in update.idx) {
        xl <- intervals[idx, xl]
        xr <- intervals[idx, xr]
        if(is.infinite(xl)) {
            # open left interval
            W <- intervals[idx, Wr]
            delta1 <- search.open.peak()
            Eu1 <- agent$slope * (delta1 + agent$peak - xr) + ct(delta1)
        } else if(is.infinite(xr)) {
            # open right interval
            W <- intervals[idx, Wl]
            delta1 <- search.open.peak()
            Eu1 <- agent$slope * (delta1 + xl - agent$peak) + ct(delta1)
        } else {
            # bridge interval
            Wl <- intervals[idx, Wl]
            Wr <- intervals[idx, Wr]
            dir <- if(xl > agent$peak) 1 else -1
            delta1 <- search.brid.peak(xl, xr, dir = dir)$maximum
            dbar <- xr - xl - delta1
            if (delta1 < 10 * .Machine$double.eps) {
                Eu1 <- -Inf
            } else {
                W1 <- if(dir > 0) { agent$slope * (delta1 + xl - agent$peak) } else { agent$slope * (delta1 - xl + agent$peak) }
                Eu1 <- W1 + ct(delta1) + ct(dbar)
            }
        }
        plan[idx, `:=`(delta=delta1, Eu=Eu1)]
    }
    agent$plan <- plan
    agent
}


