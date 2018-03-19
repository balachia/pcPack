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

############################################################
# Utility Diagnostic

#' @import ggplot2
#' @importFrom stats quantile
plot_crit <- function(critfs, n=1e3+1, ylim=0.95, xl=0, xr=1, midpoint=1, ...) {
    if(!'list' %in% class(critfs)) critfs <- list(critfs)
    xlims <- 0.5*(xr+xl + midpoint*(xr-xl)*c(-1,1))
    xs <- seq(xlims[1], xlims[2], length.out=n)
    #xs <- seq(xl, xr, length.out=n)

    ysdts <- list()
    for(i in seq_along(critfs)) {
        critf <- critfs[[i]]
        ys <- critf(xs, xl=xl, xr=xr, debug=TRUE, ...)
        cols <- c('crit', 'Mp', 'hVp', 'eurat', 'hVpeurat')
        ysdt <- rbindlist(map(cols, ~ data.table(x=xs, y=ys[[.]], type=., func=i)))
        ysdts <- c(ysdts, list(ysdt))
    }
    ysdt <- rbindlist(ysdts)

    ysdt[, maxy := quantile(abs(y), ylim, na.rm=TRUE), by=func]
    plotdt <- ysdt[abs(y) <= maxy | is.na(y) | is.infinite(y)]

    #xlims <- 0.5*(xr+xl + midpoint*(xr-xl)*c(-1,1))

    plotcols <- cols[c(2,5)]
    ggp <- ggplot(plotdt[type %in% plotcols], aes(x, y, color=type)) +
        facet_wrap(~func, ncol=1, scales='free_y') +
        geom_hline(yintercept=0) +
        #coord_cartesian(ylim=c(-ylim, ylim)) +
        coord_cartesian(expand=FALSE) +
        geom_line() +
        geom_line(data=plotdt[type=='crit'], color='black', linetype='22')
    print(ggp)

    ysdt
}

plot_bridcrit <- function(xl, xr, Wl, Wr, ...) {
    plot_crit(bridcrit, xl=xl, xr=xr, Wl=Wl, Wr=Wr, ...)
}

plot_bridcrits <- function(xl, xr, Wl, Wr, ...) {
    plot_crit(list(bridcrit.normal, bridcrit.logged),
              xl=xl, xr=xr, Wl=Wl, Wr=Wr, ...)
}

############################################################
# Search Diagnostic

# @import data.table
bridge_search_diagnostic <- function(xlseq=0, xrseq=10^seq(-5,5,length.out=1+1e2),
                                   Wlseq=0, Wrseq=10^seq(-10,5,length.out=1+1e2),
                                   ...) {
    ncol <- 50
    ngap <- 100
    grid <- CJ(xl=xlseq, xr=xrseq, Wl=Wlseq, Wr=Wrseq)
    res <- sapply(1:nrow(grid), function(i) {
        if(i %% ngap == 0) cat('.')
        if(i %% (ngap*ncol) == 0) cat(sprintf(' %3.2f%%\n', 100*i/nrow(grid)))
        tryCatch({
            xl=grid[i, xl]
            xr=grid[i, xr]
            Wl=grid[i, Wl]
            Wr=grid[i, Wr]
            res <- search.brid(xl, xr, Wl, Wr, ...)
            NA
        }, error=function(e) {
            #cat(sprintf('%s, %s: %s, %s: ', xl, xr, Wl, Wr))
            #cat('\n')
            #print(e$message)
            #print(e$call)
            e$message
        })
    })
    cat('\n')
    grid[, res := res]
    grid
}

plot_bridge_search_diagnostic <- function(dt=NULL, ...) {
    if(is.null(dt)) {
        dt <- bridge_search_diagnostic(...)
    }
    ggp <- ggplot(dt, aes(xr, Wr, fill=res)) +
        theme(legend.position='bottom') +
        scale_x_log10() +
        scale_y_log10() +
        coord_cartesian(expand=FALSE) +
        scale_fill_brewer(palette='Dark2') +
        geom_tile(alpha=0.75)
    ggp
}
