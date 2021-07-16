############################################################
# agent generics

# general flow:
# 1) create agent with `make_X_agent` and pass to simulator
# 2) `run_simulation` calls `agentEntry` on selected agent
#       expect: agent returns one entry position
# 3) `run_simulation` calls `agentUpdate` on ALL agents
#       expect: each agent updates its table of position values

#' agent entry
#'
#' Request agent's next entry position
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
#' @param agent agent to updated
#' @param update.idx updated indices in intervals table
#' @param intervals (data.table) new market intervals table
#' @param ... additional arguments to agents
#' @return updated agent
#' @export
agentUpdate <- function(agent, update.idx, intervals, ...) UseMethod('agentUpdate')

agentUpdate.default <- function(agent, update.idx, intervals, ...) agent

#' Set up agent for simulation
#'
#' @param agent agent to set up
#' @param n number of entries into market
#' @export
agentSetup <- function(agent, n, ...) UseMethod('agentSetup')

agentSetup.default <- function(agent, n, ...) agent

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
make_standard_agent <- function(a=1, b=1, ...) {
    ag <- list()
    class(ag) <- c('standard.agent', 'agent')

    ag$parameters <- list(a=a, b=b)
    ag$a <- a
    ag$b <- b
    #ag$plan <- create_plans_table(n+1)

    ag$Madj.open <- ct.Madj.open
    ag$Mpadj.open <- ct.Mpadj.open
    ag$Madj.brid <- ct.Madj.brid
    ag$Mpadj.brid <- ct.Mpadj.brid

    ag
}

#' @export
Ops.standard.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'standard.agent') && inherits(a2, 'standard.agent')) {
                   a1$a == a2$a && a1$b == a2$b
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.standard.agent <- function(ag, ...) {
    cat(sprintf('standard agent :: a %s :: b %s\n', ag$a, ag$b))
}

agentSetup.standard.agent <- function(agent, n, ...) {
    agent$plan <- create_plans_table(n+1)
    agent
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
            delta1 <- search.open(W, a=agent$a, b=agent$b,
                                  Madj=agent$Madj.open, Mpadj=agent$Mpadj.open,
                                  ...)$root
            Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b,
                          Madj=agent$Madj.open,
                          ...)
        } else if(is.infinite(xr)) {
            # open right interval
            W <- intervals[idx, Wl]
            delta1 <- search.open(W, a=agent$a, b=agent$b,
                                  Madj=agent$Madj.open, Mpadj=agent$Mpadj.open,
                                  ...)$root
            Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b,
                          Madj=agent$Madj.open,
                          ...)
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

############################################################
# GOM AGENT

#' set up entry plan table
#' 
#' Fields
#'  id --- entry sequence (corresponds to interval id)
#'  delta --- best entry in interval
#'  Eu --- expected utility of entry at delta
#' @param n number of positions
#' @import data.table
create_gom_plans_table <- function(n) {
    data.table(id=1:n, delta=as.numeric(NA), Eu=-Inf,
               cat.id=as.integer(NA), cat.mu=as.numeric(NA), cat.peak=as.numeric(NA),
               key='id')
}


#' Category-aware agent
#'
#' Agent with utility u(m) = a*m - exp(-b*m)
#' @param n number of entries into the market
#' @param a agent utility function a parameter
#' @param b agent utility function b parameter
#' @param ... capture additional arguments
#' @export
#make_gom_agent <- function(a=1, b=1, logp=TRUE, ...) {
make_gom_agent <- function(a=1, b=1, gom.weight=1, logp=FALSE, ...) {
    ag <- list()
    class(ag) <- c('gom.agent', 'standard.agent', 'agent')

    ag$parameters <- list(a=a, b=b)
    ag$a <- a
    ag$b <- b
    #ag$plan <- create_plans_table(n+1)

    ag$Madj.open <- ct.Madj.open
    ag$Mpadj.open <- ct.Mpadj.open
    ag$Madj.brid <- ct.Madj.brid
    ag$Mpadj.brid <- ct.Mpadj.brid


    # use logged goms?
    ag$logp <- logp

    # weight on gom relative to terrain
    ag$gom.weight <- gom.weight

    # build M adjustments:
    if(logp) {
        ag$gom.Madj <- function(...) gom.weight*gom.log.Madj(...)
        ag$gom.Mpadj <- function(...) gom.weight*gom.log.Mpadj(...)
    } else {
        ag$gom.Madj <- function(...) gom.weight*gom.Madj(...)
        ag$gom.Mpadj <- function(...) gom.weight*gom.Mpadj(...)
    }

    ag
}

#' @export
Ops.gom.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'gom.agent') && inherits(a2, 'gom.agent')) {
                   a1$a == a2$a && a1$b == a2$b
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.gom.agent <- function(ag, ...) {
    logflag <- if(ag$logp) { 'log' } else { 'nolog' }
    cat(sprintf('gom agent :: a %s :: b %s :: gw %s :: %s\n',
                ag$a, ag$b, ag$gom.weight, logflag))
}

agentSetup.gom.agent <- function(agent, n, ...) {
    agent$plan <- create_gom_plans_table(n+1)
    agent
}

#' @import data.table
agentUpdate.gom.agent <- function(agent, update.idx, intervals, verbose=.verbose$NONE, ...) {
    extra <- list(...)
    debug <- if(!is.null(extra$debug)) { extra$debug } else { FALSE }
    plan <- agent$plan
    #gw <- agent$gom.weight
    # rebuild categories
    xlfin <- intervals[, is.finite(xl)]
    xrfin <- intervals[, is.finite(xr)]
    either <- xlfin|xrfin
    xs <- intervals[(xlfin), xl]
    # only do this categorization if we have enough data to categorize
    if(length(xs) >= 3) {
        #   get mixture
        extra$verbose <- FALSE
        mix <- do.call(categorize.mclust, c(list(xs, min.k=2), extra))$mix
        #mix <- categorize.mclust(xs, min.k=2, verbose=FALSE, ...)$mix
        mixps <- get_mix_parameters(mix)
        peaks <- sapply(1:mixps$k, function(i) {
                dnorm(mixps$mean[i], mean=mixps$mean[i], sd=mixps$sd[i], log=TRUE)
            })
        #   get goms & top category
        xgoms <- intervals[(either), ifelse(xlfin[either], xl, xr)]
        goms <- get_goms(xgoms, mix, logp=TRUE)
        topcat <- sapply(1:nrow(goms), function(ri) order(goms[ri,], decreasing=TRUE)[1])
        plan[(either), `:=`(cat.id=topcat,
                                 cat.mu=mixps$mean[topcat],
                                 cat.peak=peaks[topcat])]

        # find valid intervals
        idxs <- intervals[(either), id]
        # predraw vectors
        xls <- intervals[, xl]
        xrs <- intervals[, xr]
        Wls <- intervals[, Wl]
        Wrs <- intervals[, Wr]
        cat.ids <- plan[, cat.id]
        # preallocate results
        deltas <- numeric(length=length(idxs))
        Eus <- numeric(length=length(idxs))
        # update each interval
        for(idx in idxs) {
            #xl <- intervals[idx, xl]
            #xr <- intervals[idx, xr]
            #idx.cat <- plan[idx, cat.id]
            xl <- xls[idx]
            xr <- xrs[idx]
            idx.cat <- cat.ids[idx]
            gom.mean <- mixps$mean[idx.cat]
            gom.sd <- mixps$sd[idx.cat]
            gom.peak <- peaks[idx.cat]
            if(debug) {
                cat(sprintf('gom looking at %d :: x (%s, %s) W (%s, %s), mu %s sd %s\n',
                                  idx, xl, xr, intervals[idx, Wl], intervals[idx, Wr],
                                  gom.mean, gom.sd))
            }
            if(is.infinite(xl)) {
                # open left interval
                Madj <- function(delta, ...) agent$Madj.open(delta) + agent$gom.Madj(xr-delta, gom.mean, gom.sd, gom.peak, ...)
                Mpadj <- function(delta, ...) agent$Mpadj.open(delta) - agent$gom.Mpadj(xr-delta, gom.mean, gom.sd, gom.peak, ...)
                #W <- intervals[idx, Wr]
                W <- Wrs[idx]
                delta1 <- search.open(W, a=agent$a, b=agent$b,
                                      Madj=Madj, Mpadj=Mpadj, verbose=verbose,
                                      ...)$root
                Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b,
                              Madj=Madj, verbose=verbose,
                              ...)
            } else if(is.infinite(xr)) {
                # open right interval
                Madj <- function(delta, ...) agent$Madj.open(delta) + agent$gom.Madj(xl+delta, gom.mean, gom.sd, gom.peak, ...)
                Mpadj <- function(delta, ...) agent$Mpadj.open(delta) + agent$gom.Mpadj(xl+delta, gom.mean, gom.sd, gom.peak, ...)
                #W <- intervals[idx, Wl]
                W <- Wls[idx]
                delta1 <- search.open(W, a=agent$a, b=agent$b,
                                      Madj=Madj, Mpadj=Mpadj, verbose=verbose,
                                      ...)$root
                Eu1 <- Euopen(delta1, W0=W, a=agent$a, b=agent$b,
                              Madj=Madj, verbose=verbose,
                              ...)
            } else {
                # bridge interval
                Madj <- function(delta, dbar, ...) agent$Madj.brid(delta, dbar) + agent$gom.Madj(xl+delta, gom.mean, gom.sd, gom.peak, ...)
                Mpadj <- function(delta, dbar, ...) agent$Mpadj.brid(delta, dbar) + agent$gom.Mpadj(xl+delta, gom.mean, gom.sd, gom.peak, ...)
                #Wl <- intervals[idx, Wl]
                #Wr <- intervals[idx, Wr]
                Wl <- Wls[idx]
                Wr <- Wrs[idx]
                delta1 <- search.brid(xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b,
                                      Madj=Madj, Mpadj=Mpadj, verbose=verbose,
                                      ...)$u0
                Eu1 <- Eubrid(delta1, xl=xl, xr=xr, Wl=Wl, Wr=Wr, a=agent$a, b=agent$b,
                              Madj=Madj, verbose=verbose,
                              ...)
            }
            #plan[idx, `:=`(delta=delta1, Eu=Eu1)]
            deltas[idx] <- delta1           
            Eus[idx] <- Eu1           
        }
        plan[idxs, `:=`(delta=deltas, Eu=Eus)]
        agent$plan <- plan
    }
    agent
}

############################################################
# INSERT AGENT

make_insert_agent <- function(insert.dt, ...) {
    ag <- list()
    class(ag) <- c('insert.agent', 'agent')

    ag$parameters <- list(insert.dt=insert.dt)
    ag$insert.dt <- insert.dt

    ag
}

#' @export
print.insert.agent <- function(ag, ...) {
    cat(sprintf('insert agent :: n\n', nrow(ag$insert.dt)))
}

#' @import data.table
agentEntry.insert.agent <- function(agent, i, intervals, positions, ...) {
    x <- agent$insert.dt[i, x]
    W <- agent$insert.dt[i, W]

    enteridx <- intervals[x > xl & x < xr, id[1]]
    list(idx=enteridx, x=x, W=W)
}
