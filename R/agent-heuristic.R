############################################################
# heuristic agent

#' Standard utility function agent
#'
#' Heuristic agent that always enters at midpoint of bridges and at `jump` distance of open intervals
#' @param n number of entries into the market
#' @param jump agent jump distance on open intervals
#' @param sigma landscape sigma parameter
#' @param ... capture additional arguments
#' @export
make_heuristic_agent <- function(jump=1, sigma = 1, ...) {
    ag <- list()
    class(ag) <- c('heuristic.agent', 'agent')

    ag$parameters <- list(jump = jump, sigma = sigma)
    ag$jump <- jump
    ag$sigma <- sigma
    #ag$plan <- create_plans_table(n+1)

    ag$Madj.open <- ct.Madj.open
    ag$Mpadj.open <- ct.Mpadj.open
    ag$Madj.brid <- ct.Madj.brid
    ag$Mpadj.brid <- ct.Mpadj.brid

    ag
}

#' @export
Ops.heuristic.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'heuristic.agent') && inherits(a2, 'heuristic.agent')) {
                   a1$jump == a2$jump && a1$sigma == a2$sigma
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.heuristic.agent <- function(ag, ...) {
    cat(sprintf('heuristic agent :: jump %s :: sd %s\n', ag$jump, ag$sigma))
}

agentSetup.heuristic.agent <- function(agent, n, ...) {
    agent$plan <- create_plans_table(n+1)
    agent
}

#' @import data.table
agentEntry.heuristic.agent <- function(agent, i, intervals, positions, ...) {
    enteridx <- agent$plan[which.max(Eu), id]
    delta <- agent$plan[enteridx, delta]

    # entering at open or closed interval?
    xl <- intervals[enteridx, xl]
    xr <- intervals[enteridx, xr]
    if(is.infinite(xl) || is.infinite(xr)) {
        Wd <- openjump(delta, sigma = agent$sigma, ...)
        if(is.infinite(xl)) {
            x <- xr-delta
            W <- intervals[enteridx, Wr+Wd]
        } else {
            x <- xl+delta
            W <- intervals[enteridx, Wl+Wd]
        }
    } else {
        x <- xl+delta
        W <- bridjump(delta, xl=xl, xr=xr,
                      sigma = agent$sigma,
                      Wl=intervals[enteridx, Wl],
                      Wr=intervals[enteridx, Wr], ...)
    }

    list(idx=enteridx, x=x, W=W)
}

#' @import data.table
agentUpdate.heuristic.agent <- function(agent, update.idx, intervals, ...) {
    plan <- agent$plan
    for(idx in update.idx) {
        xl <- intervals[idx, xl]
        xr <- intervals[idx, xr]
        if(is.infinite(xl)) {
            # open left interval
            W <- intervals[idx, Wr]
            delta1 <- agent$jump
            Eu1 <- W + ct(delta1)
        } else if(is.infinite(xr)) {
            # open right interval
            W <- intervals[idx, Wl]
            delta1 <- agent$jump
            Eu1 <- W + ct(delta1)
        } else {
            # bridge interval
            Wl <- intervals[idx, Wl]
            Wr <- intervals[idx, Wr]
            delta1 <- (xr - xl) / 2
            if (delta1 < 10 * .Machine$double.eps) {
                Eu1 <- -Inf
            } else {
                Eu1 <- ((Wl + Wr) / 2) + (2 * ct(delta1))
            }
        }
        plan[idx, `:=`(delta=delta1, Eu=Eu1)]
    }
    agent$plan <- plan
    agent
}


