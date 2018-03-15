############################################################
# Simulation Diagnostic and Profiling

#' @import utils
profile_simulation <- function(n=1e2, vis=TRUE) {
    tmp <- tempfile()

    set.seed(1)

    Rprof(tmp, line.profiling=TRUE)
    test_simulation(n)
    Rprof(NULL)

    summaryRprof(tmp, lines='show')

    if(vis) {
        profvis::profvis(prof_input=tmp)
    }
}

profvis_simulation <- function(n=1e2) {
    profvis::profvis({
        test_simulation(n)
    })
}

test_simulation <- function(n=100, diagnose.agent=FALSE) {
    insert.dt <- data.table(x=0, W=0)
    if(diagnose.agent) {
        ag <- make_diagnostic_agent(n)
    } else {
        ag <- make_standard_agent(n)
    }

    agl <- set_up_agents(n, insert.dt,
                         list(ag),
                         randomize=FALSE)

    run_simulation(n, agl$agents, agl$order)
}

############################################################
# Agent Diagnostic

make_diagnostic_agent <- function(...) {
    ag <- make_standard_agent(...)
    class(ag) <- c('diagnostic.agent', class(ag))
    ag
}

#' @import ggplot2
agentEntry.diagnostic.agent <- function(ag, i, intervals, positions, ...) {
    plan <- ag$plan[, list(id, delta, Eu=-log1p(max(Eu) - Eu))]
    res <- agentEntry.standard.agent(ag, i, intervals, positions, ...)

    #plotdt <- intervals[, list(id, x=xl)]
    plotdt <- positions[!is.na(x), list(x, y=W, type='W')]
    plotdt <- rbind(plotdt,
                    intervals[plan[is.finite(Eu)],
                              list(x=ifelse(is.finite(xl), xl+delta, xr-delta),
                                   y=Eu,
                                   type='Eu')])

    ggp <- ggplot(plotdt, aes(x, y)) +
        facet_wrap(~type, ncol=1, scales='free_y') +
        geom_vline(xintercept=res$x) +
        geom_point() +
        geom_point(x=res$x, y=res$W, color='red')
    print(ggp)

    #readline()
    #Sys.sleep(0.1)

    res
}

############################################################
# Brownian Walk Diagnostic

#' @import ggplot2
#' @importFrom utils head tail
diagnose_positions <- function(positions) {
    dt <- positions[order(x),
                    list(dx = tail(x, -1) - head(x, -1),
                         dW = tail(W, -1) - head(W, -1))]
    dt[, z := dW/sqrt(dx)]

    ggp <- ggplot(dt, aes(x=z)) +
        geom_histogram(bins=10)
    #print(ggp)

    ggp <- ggplot(dt, aes(sample=z)) +
        coord_fixed() +
        geom_abline(slope=1, intercept=0) +
        stat_qq()
    print(ggp)

    dt
}
