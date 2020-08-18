############################################################
##### Salop circle agent
# Assumes operating on a 'circular' brownian bridge, bounded by outermost positions.
# Will never jump to open intervals.
# Mechanically: standard agent, but treats open intervals as Eu = -Inf ??

#' Standard utility function agent on Salop circle model
#'
#' Agent with utility u(m) = a*m - exp(-b*m)
#'
#' @param n number of entries into the market
#' @param a agent utility function a parameter
#' @param b agent utility function b parameter
#' @param ... capture additional arguments
#' @export
make_salop_agent <- function(a=1, b=1, ...) {
    ag <- make_standard_agent(a=a, b=b, ...)
    #class(ag) <- c('salop.agent', 'standard.agent', 'agent')
    class(ag) <- c('salop.agent', 'agent')

    ag
}

#' @export
Ops.salop.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'salop.agent') && inherits(a2, 'salop.agent')) {
                   a1$a == a2$a && a1$b == a2$b
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.salop.agent <- function(ag, ...) {
    cat(sprintf('salop standard agent :: a %s :: b %s\n', ag$a, ag$b))
}

agentSetup.salop.agent <- function(agent, n, ...) {
    agent$plan <- create_plans_table(n+1)
    agent
}

#' @import data.table
agentEntry.salop.agent <- function(agent, i, intervals, positions, ...) {
    agentEntry.standard.agent(agent, i, intervals, positions, ...)
    # return: list(idx=enteridx, x=x, W=W)
}

#' @import data.table
agentUpdate.salop.agent <- function(agent, update.idx, intervals, ...) {
    plan <- agent$plan
    for(idx in update.idx) {
        xl <- intervals[idx, xl]
        xr <- intervals[idx, xr]
        if(is.infinite(xl)) {
            # open left interval: salop never accepts open interval
            W <- intervals[idx, Wr]
            delta1 <- 0
            Eu1 <- -Inf
        } else if(is.infinite(xr)) {
            # open right interval: salop never accepts open interval
            delta1 <- 0
            Eu1 <- -Inf
        } else {
            # bridge interval
            Wl <- intervals[idx, Wl]
            Wr <- intervals[idx, Wr]
            delta1 <- search.brid(xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b,
                                  Madj=agent$Madj.brid, Mpadj=agent$Mpadj.brid,
                                  ...)$u0
            Eu1 <- Eubrid(delta1, xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b,
                          Madj=agent$Madj.brid,
                          ...)
        }
        plan[idx, `:=`(delta=delta1, Eu=Eu1)]
    }
    agent$plan <- plan
    agent
}



