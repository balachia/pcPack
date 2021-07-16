############################################################
##### Agent operating under softmax/logistic competition rules
# brownian bridge represents 'quality' of good, competing against other producers and reserve good at q = 0
# consumer at x selects good at x0, valued at q0, proportional to exp(q0 - |x - x0|)
# this leads to a market share, and producers optimize their positions subject
#   to a utility function that accounts for the option value of selling trash
#   under high variance exploration

# TODO:
#' Softmax utility agent
#'
#' Agent with utility u(m) = log(exp(a*(m/2)) - 1) -- a certain inverse of logistic market share
#' @param n number of entries into the market
#' @param a agent utility function a parameter
#' @param ... capture additional arguments
#' @export
make_softmax_agent <- function(a=1, ...) {
    ag <- list()
    class(ag) <- c('softmax.agent', 'agent')

    ag$parameters <- list(a=a)
    ag$a <- a
    ag$plan <- list(x0 = NA_real_,
                    enteridx = NA_integer_,
                    delta = NA_real_,
                    enter.joint = NA,
                    term.optim = c(-4, 4))

    ag
}

# TODO:
#' @export
Ops.softmax.agent <- function(a1, a2) {
    switch(.Generic[[1]],
           `==` = {
               if(inherits(a1, 'softmax.agent') && inherits(a2, 'softmax.agent')) {
                   a1$a == a2$a
               } else {
                   FALSE
               }
           },
           stop('Agent operator not implemented')
           )
}

#' @export
print.softmax.agent <- function(ag, ...) {
    cat(sprintf('softmax agent :: a %s\n', ag$a))
}

# TODO:
#' @importFrom fastGHQuad gaussHermiteData
agentSetup.softmax.agent <- function(agent, gh.n = 30, ...) {
    #agent$plan <- create_plans_table(n+1)
    agent$gh.rule <- fastGHQuad::gaussHermiteData(gh.n)
    agent
}

# TODO:
#' @import data.table
#' @importFrom pcSoftmaxPack Eumsi.optimize
agentEntry.softmax.agent <- function(agent, i, intervals, positions, ...) {
    # get agent plan
    x <- agent$plan$x0
    enteridx <- agent$plan$enteridx
    delta <- agent$plan$delta

    #enteridx <- agent$plan[which.max(Eu), id]
    #delta <- agent$plan[enteridx, delta]

    # entering at open or closed interval?
    xl <- intervals[enteridx, xl]
    xr <- intervals[enteridx, xr]
    if(is.infinite(xl) || is.infinite(xr)) {
        Wd <- openjump(delta, ...)
        if(is.infinite(xl)) {
            #x <- xr-delta
            W <- intervals[enteridx, Wr+Wd]
        } else {
            #x <- xl+delta
            W <- intervals[enteridx, Wl+Wd]
        }
    } else {
        #x <- xl+delta
        W <- bridjump(delta, xl=xl, xr=xr,
                      Wl=intervals[enteridx, Wl],
                      Wr=intervals[enteridx, Wr], ...)
    }

    list(idx=enteridx, x=x, W=W)
}

#' @import data.table
agentUpdate.softmax.agent <- function(agent, update.idx, intervals, ...) {
    # this agent does not need to make plans:
    # everything relevant can be calculated at entry and every interval needs to be recalculated anyway
    # at least as far as i can tell
    # it's probably possible to sort 'likely' best locations, but not sure what the guarantees are

    # check debug
    debug <- getOption('debug', FALSE)

    # reconstruct xgrid, vgrid from intervals
    # make sure to remember associations for which interval ppl enter at
    # find valid interval endpoints
    xlfin <- intervals[, is.finite(xl)]
    xrfin <- intervals[, is.finite(xr)]
    terminus <- xor(xlfin, xrfin)
    either <- xlfin|xrfin
    # grid consists of non-infinite left endpoints (equivalent to non-infinite right endpoints)
    xgrid <- intervals[(xlfin), xl]
    vgrid <- intervals[(xlfin), Wl]

    #if(debug) cat(sprintf('softmax agent starting plan\n'))

    # run the optimization
    # returns optimum position, and whether it is terminal, joint, or otherwise
    optima <- pcSoftmaxPack::Eumsi.optimize(vgrid = vgrid, xgrid = xgrid,
                                            term.last.optim = agent$plan$term.optim,
                                            gh.rule = agent$gh.rule, a = agent$a)

    #if(debug) cat(sprintf('softmax agent done\n'))

    x0 <- optima$maximum
    enter.joint <- optima$is.joint
    if(optima$is.joint) {
        # if optimum is at a joint, insert at the unique interval which contains it as a left but not right endpoint
        # find that interval, but make sure we're not hitting the infinite intervals
        enteridx <- intervals[either & xl <= x0 & xr > x0, id]
        if(length(enteridx) != 1) {
            print(optima)
            print(enteridx)
            stop("Softmax agent failed at entering joint")
        }
        stopifnot(length(enteridx) == 1)
        #enteridx <- jointidx
        delta <- 0
    } else if(optima$is.terminal) {
        # try left terminus
        eidxl <- intervals[!(xlfin) & (xrfin) & (x0 <= xr), .(idx = id, delta = xr - x0)]
        # try right terminus
        eidxr <- intervals[(xlfin) & !(xrfin) & (x0 >= xl), .(idx = id, delta = x0 - xl)]
        enteridx <- c(eidxl$idx, eidxr$idx)
        if(length(enteridx) != 1) {
            print(optima)
            print(enteridx)
            stop("Softmax agent failed at entering terminus")
        }
        #enteridx <- termidx
        delta <- c(eidxl$delta, c(eidxr$delta))
    } else {
        enteridx <- intervals[either & xl <= x0 & xr > x0, id]
        if(length(enteridx) != 1) {
            print(optima)
            print(enteridx)
            stop("Softmax agent failed at entering interval")
        }
        #enteridx <- itlidx
        delta <- intervals[enteridx, x0 - xl]
    }

    if(debug) {
        cat(sprintf('softmax agent plans: %f (itl %d), jumping %f\n', x0, enteridx, delta))
    }

    agent$plan <- list(x0 = x0, enteridx = enteridx, delta = delta,
                       enter.joint = enter.joint, term.optim = optima$term.optim)
    agent
}
