############################################################
# agent generics

#' agent entry
#'
#' @param agent entering agent
#' @param i entry count (id)
#' @param intervals (data.table) existing market intervals table
#' @param positions (data.table) existing market positions table
#' @param ... additional arguments to agents
#' @export
agentEntry <- function(agent, i, intervals, positions, ...) UseMethod('agentEntry')

#' agent update
#'
#' @param agent updated agent
#' @param update.idx updated indices in intervals table
#' @param intervals (data.table) new market intervals table
#' @param ... additional arguments to agents
#' @export
agentUpdate <- function(agent, update.idx, intervals, ...) UseMethod('agentUpdate')

agentUpdate.default <- function(agent, update.idx, intervals, ...) agent

############################################################
# STANDARD AGENT

#' set up entry plan table
#' 
#' Fields
#'  id --- entry sequence (corresponds to interval id)
#'  delta --- best entry in interval
#'  Eu --- expected utility of entry at delta
#' @param n number of positions
#' @import data.table
create_plans_table <- function(n) {
    data.table(id=1:n, delta=as.numeric(NA), Eu=-Inf, key='id')
}

#' Standard utility function agent
#'
#' Agent with utility u(m) = a*m - exp(-b*m)
#' @param n number of entries into the market
#' @param a agent utility function a parameter
#' @param b agent utility function b parameter
#' @param ... capture additional arguments
#' @export
make_standard_agent <- function(n, a=1, b=1, ...) {
    ag <- list()
    class(ag) <- c('standard.agent', 'agent')

    ag$parameters <- list(a=a, b=b)
    ag$a <- a
    ag$b <- b
    ag$plan <- create_plans_table(n+1)

    ag
}

#' @import data.table
agentEntry.standard.agent <- function(agent, i, intervals, positions, ...) {
    enteridx <- agent$plan[which.max(Eu), id]
    delta <- agent$plan[enteridx, delta]

    # entering at open or closed interval?
    xl <- intervals[enteridx, xl]
    xr <- intervals[enteridx, xr]
    if(is.infinite(xl) || is.infinite(xr)) {
        Wd <- openjump(delta, ...)
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
                      Wl=intervals[enteridx, Wl],
                      Wr=intervals[enteridx, Wr], ...)
    }

    list(idx=enteridx, x=x, W=W)
}

#' @import data.table
agentUpdate.standard.agent <- function(agent, update.idx, intervals, ...) {
    plan <- agent$plan
    for(idx in update.idx) {
        xl <- intervals[idx, xl]
        xr <- intervals[idx, xr]
        if(is.infinite(xl)) {
            # open left interval
            W <- intervals[idx, Wr]
            delta1 <- search.open(W, a=agent$a, b=agent$b, ...)$root
            Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b, ...)
        } else if(is.infinite(xr)) {
            # open right interval
            W <- intervals[idx, Wl]
            delta1 <- search.open(W, a=agent$a, b=agent$b, ...)$root
            Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b, ...)
        } else {
            # bridge interval
            Wl <- intervals[idx, Wl]
            Wr <- intervals[idx, Wr]
            delta1 <- search.brid(xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b, ...)$u0
            Eu1 <- Eubrid(delta1, xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b, ...)
        }
        plan[idx, `:=`(delta=delta1, Eu=Eu1)]
    }
    agent$plan <- plan
    agent
}

############################################################
# INSERT AGENT

make_insert_agent <- function(n, insert.dt, ...) {
    ag <- list()
    class(ag) <- c('insert.agent', 'agent')

    ag$parameters <- list(insert.dt=insert.dt)
    ag$insert.dt <- insert.dt

    ag
}

#' @import data.table
agentEntry.insert.agent <- function(agent, i, intervals, positions, ...) {
    x <- agent$insert.dt[i, x]
    W <- agent$insert.dt[i, W]

    enteridx <- intervals[x > xl & x < xr, id[1]]
    list(idx=enteridx, x=x, W=W)
}
