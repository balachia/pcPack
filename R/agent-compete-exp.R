############################################################
# STANDARD AGENT, competition function exponential in competitor value

############################################################
# COMPETITION 

exp.compete.competition.functions <- function() {
    #' Competition function
    #'
    #' Competition function
    #' @param delta distance from competitor
    ct.exp <- function(delta, W, ...) -exp(W) / delta

    #' Competition function (d/ddelta)
    #'
    #' Competition function (d/ddelta)
    #' @param delta distance from competitor
    dct.exp <- function(delta, W, ...) exp(W) / (delta^2)

    Madj.open <- function(delta, W, ...) ct.exp(delta, W)
    Mpadj.open <- function(delta, W, ...) dct.exp(delta, W)
    Madj.brid <- function(delta, dbar, Wl, Wr, ...) ct.exp(delta, W = Wl) + ct.exp(dbar, W = Wr)
    Mpadj.brid <- function(delta, dbar, Wl, Wr, ...) dct.exp(delta, W = Wl) - dct.exp(dbar, W = Wr)
    list(Madj.open = Madj.open, Mpadj.open = Mpadj.open, Madj.brid = Madj.brid, Mpadj.brid = Mpadj.brid)
}


############################################################
# agent definition

#' Standard utility function agent
#'
#' Agent with utility u(m) = a*m - exp(-b*m)
#' @param n number of entries into the market
#' @param a agent utility function a parameter
#' @param b agent utility function b parameter
#' @param ... capture additional arguments
#' @export
make_exp_compete_agent <- function(a=1, b=1, sigma = 1, ...) {
    ag <- list()
    class(ag) <- c('standard.exp.compete.agent', 'agent')

    ag$parameters <- list(a=a, b=b)
    ag$a <- a
    ag$b <- b
    ag$sigma <- sigma
    #ag$plan <- create_plans_table(n+1)

    # TODO: change these
    cts <- exp.compete.competition.functions()
    ag$Madj.open <- cts$Madj.open
    ag$Mpadj.open <- cts$Mpadj.open
    ag$Madj.brid <- cts$Madj.brid
    ag$Mpadj.brid <- cts$Mpadj.brid

    ag
}

#' @export
Ops.standard.exp.compete.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'standard.exp.compete.agent') && inherits(a2, 'standard.exp.compete.agent')) {
                   a1$a == a2$a && a1$b == a2$b
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.standard.exp.compete.agent <- function(ag, ...) {
    cat(sprintf('standard exp-compete agent :: a %s :: b %s :: sd %s\n', ag$a, ag$b, ag$sigma))
}

agentSetup.standard.exp.compete.agent <- function(agent, n, ...) {
    agent$plan <- create_plans_table(n+1)
    agent
}

#' @import data.table
agentEntry.standard.exp.compete.agent <- function(agent, i, intervals, positions, ...) {
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
agentUpdate.standard.exp.compete.agent <- function(agent, update.idx, intervals, ...) {
    plan <- agent$plan
    for(idx in update.idx) {
        xl <- intervals[idx, xl]
        xr <- intervals[idx, xr]
        if(is.infinite(xl)) {
            # open left interval
            W <- intervals[idx, Wr]
            delta1 <- search.open(W, a=agent$a, b=agent$b,
                                  Madj=agent$Madj.open, Mpadj=agent$Mpadj.open,
                                  sigma = agent$sigma,
                                  ...)$root
            Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b,
                          Madj=agent$Madj.open,
                          sigma = agent$sigma,
                          ...)
        } else if(is.infinite(xr)) {
            # open right interval
            W <- intervals[idx, Wl]
            delta1 <- search.open(W, a=agent$a, b=agent$b,
                                  Madj=agent$Madj.open, Mpadj=agent$Mpadj.open,
                                  sigma = agent$sigma,
                                  ...)$root
            Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b,
                          Madj=agent$Madj.open,
                          sigma = agent$sigma,
                          ...)
        } else {
            # bridge interval
            Wl <- intervals[idx, Wl]
            Wr <- intervals[idx, Wr]
            delta1 <- search.brid(xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b,
                                  Madj=agent$Madj.brid, Mpadj=agent$Mpadj.brid,
                                  sigma = agent$sigma,
                                  ...)$u0
            Eu1 <- Eubrid(delta1, xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b,
                          Madj=agent$Madj.brid,
                          sigma = agent$sigma,
                          ...)
        }
        plan[idx, `:=`(delta=delta1, Eu=Eu1)]
    }
    agent$plan <- plan
    agent
}


