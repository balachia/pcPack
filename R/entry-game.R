#' Create market information table
#'
#' market consists of a sequence of enterable intervals
#' Fields
#'  id --- (key) entry sequence
#'  xl, xr --- location of left (right) end; +/-Inf for open intervals
#'  Wl, Wr --- value at left and right intervals; NA for unentered and open intervals
#'  state --- interval state
#' @param n number of intervals to create (n.b. should 1+number of entries)
create.intervals.table <- function(n) {
    data.table(id=1:n, xl=-Inf, xr=Inf,
               Wl=as.numeric(NA), Wr=as.numeric(NA),
               key='id')
}

#' Create market entries table
#'
#' Fields
#'  id --- (key) entry sequence
#'  x --- entry position
#'  W --- position value
#' @param n number of positions to create
create.positions.table <- function(n) {
    data.table(id=1:n, x=as.numeric(NA), W=as.numeric(NA), agent.id=as.numeric(NA),
               key='id')
}

#' Run simulation with specified agents
#'
#' @param n number of entries into the market (simulation iterations)
#' @param agents list of agent types in the market
#' @param agent.order agent entry order into the market
#' @param verbose verbosity
#' @param ... additional arguments to agents
#' @import data.table
#' @export
run_simulation <- function(n, agents,
                       agent.order=rep(seq_along(agents), length.out=n),
                       verbose=.verbose$NONE,
                       ...) {

    # prepare tables
    # intervals table is one larger than number of intervals
    # since we start with a greenfield interval open on both sides
    positions <- create.positions.table(n)
    intervals <- create.intervals.table(n+1)

    verb.interval <- n/20

    for(i in 1:n) {
        #print(intervals)
        if (verbose >= .verbose$TRACE) {
            cat(sprintf('Market iteration: %d\n', i)) 
        } else if(verbose >= .verbose$INFO && i %% verb.interval < 1) {
            cat('.')
        }
        #if(verbose >= .verbose$DEBUG && i %% verb.interval < 1) cat('.')

        agent <- agents[[agent.order[i]]]
        entry <- agentEntry(agent, i, intervals, positions, verbose=verbose, ...)

        #cat(sprintf('%s --- %s --- %s\n', entry$idx, entry$x, entry$W))

        xe <- entry$x
        We <- entry$W
        enteridx <- entry$idx

        # update tables
        positions[i, `:=`(x=xe, W=We, agent.id=agent.order[i])]
        intervals[i+1, `:=`(xl=xe,
                            xr=intervals[enteridx, xr],
                            Wl=We,
                            Wr=intervals[enteridx, Wr])]
        intervals[enteridx, `:=`(xr=xe, Wr=We)]

        # update all agents
        for(agent.i in seq_along(agents)) {
            agent <- agents[[agent.i]]
            agents[[agent.i]] <- agentUpdate(agent, c(enteridx, i+1), intervals, verbose=verbose, ...)
        }
    }

    list(positions=positions, intervals=intervals, agents=agents)
}

#' Run multiple simulations
#' 
#' @param nsim number of markets to simulation
#' @param n number of entries into each market (simulation iterations)
#' @param agent.fs list of functions to create agent types
#' @param insert.dt list of initial positions to insert
#' @param seed control random seed
#' @param verbose verbosity
#' @param ... additional arguments for agent types
#' @import data.table
#' @export
run_simulations <- function(nsim, n,
                            #agent.fs=list(make_standard_agent),
                            agents=list(make_standard_agent()),
                            insert.dt=data.table(x=0, W=0),
                            seed=1, verbose=.verbose$NONE,
                            ...) {
    # pre-sample seeds for parallel computation
    set.seed(seed)
    seeds <- sample.int(.Machine$integer.max, nsim)
    sim.format <- sprintf('SIM %%%dd / %d ', 1+floor(log10(nsim)), nsim)
    time.format <- ' (%0.2fs | %0.2fs/i)\n'
    ptm <- proc.time()
    ress <- parallel::mclapply(1:nsim,
        FUN=function(simi) {
            set.seed(seeds[simi])
            if(verbose >= .verbose$NONE) cat(sprintf(sim.format, simi))
            if(verbose >= .verbose$INFO) cat(sprintf('(seed %s) ', seeds[simi]))
            #ags0 <- lapply(agent.fs, function(f) f(n, ...))
            #agl <- set_up_agents(n, insert.dt, ags0, ...)
            agl <- set_up_agents(n, insert.dt, agents, ...)
            res <- run_simulation(n, agl$agents, agl$order, verbose=verbose, ...)
            dtime <- (proc.time() - ptm)[3]
            if(verbose >= .verbose$NONE) cat(sprintf(time.format, dtime, dtime/simi))
            res$positions[, sim := simi]
            res$intervals[, sim := simi]
            res
        })

    purrr::transpose(ress)
}

#' Combine list of simulations into a single structure
#' 
#' @param sims_list list of simulations
#' @export
combine_simulations <- function(sims_list) {
    positions <- rbindlist(sims_list$positions)
    intervals <- rbindlist(sims_list$intervals)
    list(positions=positions, intervals=intervals, agents=sims_list$agents)
}

set_up_agents <- function(n, insert.dt, agents, randomize=FALSE, ...) {
    nagents <- length(agents)
    ninsert <- nrow(insert.dt)
    #agents1 <- c(list(make_insert_agent(n, insert.dt)), agents)
    inserter <- make_insert_agent(insert.dt)
    agents1 <- lapply(c(list(inserter), agents), function(ag) agentSetup(ag, n, ...))
    if(randomize) {
        agent.order <- sample(1+(1:nagents), size=n-ninsert, replace=TRUE)
    } else {
        agent.order <- rep(1+(1:nagents), length.out=n-ninsert)
    }
    agent.order <- c(rep(1, ninsert), agent.order)

    list(agents=agents1, order=agent.order)
}
