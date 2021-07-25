############################################################
# STANDARD AGENT, variable sigma

#' Standard utility function agent
#'
#' Agent with utility u(m) = a*m - exp(-b*m)
#' @param n number of entries into the market
#' @param a agent utility function a parameter
#' @param b agent utility function b parameter
#' @param ... capture additional arguments
#' @export
make_standard_sigma_agent <- function(a=1, b=1, sigma = 1, ...) {
    ag <- list()
    class(ag) <- c('standard.sigma.agent', 'agent')

    ag$parameters <- list(a=a, b=b)
    ag$a <- a
    ag$b <- b
    ag$sigma <- sigma
    #ag$plan <- create_plans_table(n+1)

    ag$Madj.open <- ct.Madj.open
    ag$Mpadj.open <- ct.Mpadj.open
    ag$Madj.brid <- ct.Madj.brid
    ag$Mpadj.brid <- ct.Mpadj.brid

    ag
}

#' @export
Ops.standard.sigma.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'standard.sigma.agent') && inherits(a2, 'standard.sigma.agent')) {
                   a1$a == a2$a && a1$b == a2$b
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.standard.sigma.agent <- function(ag, ...) {
    cat(sprintf('standard sigma agent :: a %s :: b %s :: sd %s\n', ag$a, ag$b, ag$sigma))
}

agentSetup.standard.sigma.agent <- function(agent, n, ...) {
    agent$plan <- create_plans_table(n+1)
    agent
}

#' @import data.table
agentEntry.standard.sigma.agent <- function(agent, i, intervals, positions, ...) {
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
agentUpdate.standard.sigma.agent <- function(agent, update.idx, intervals, ...) {
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

