#' Create market information table
#'
#' market consists of a sequence of enterable intervals
#' Fields
#'  id --- (key) entry sequence
#'  xl, xr --- location of left (right) end; +/-Inf for open intervals
#'  Wl, Wr --- value at left and right intervals; NA for unentered and open intervals
#'  state --- interval state
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
create.positions.table <- function(n) {
    data.table(id=1:n, x=as.numeric(NA), W=as.numeric(NA), agent.id=as.numeric(NA),
               key='id')
}

#' Run simulation with specified agents
#'
#' @import data.table
#' @export
run_simulation <- function(n, agents,
                       agent.order=rep(seq_along(agents), length.out=n),
                       ...) {

    # prepare tables
    # intervals table is one larger than number of intervals
    # since we start with a greenfield interval open on both sides
    positions <- create.positions.table(n)
    intervals <- create.intervals.table(n+1)

    for(i in 1:n) {
        #print(intervals)

        agent <- agents[[agent.order[i]]]
        entry <- agentEntry(agent, i, intervals, positions, ...)

        cat(sprintf('%s --- %s --- %s\n', entry$idx, entry$x, entry$W))

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
            agents[[agent.i]] <- agentUpdate(agent, c(enteridx, i+1), intervals)
        }
    }

    list(positions=positions, intervals=intervals, agents=agents)
}

set_up_agents <- function(n, insert.dt, agents, randomize=FALSE) {
    agents1 <- c(list(make_insert_agent(n, insert.dt)), agents)
    if(randomize) {
        agent.order <- sample(1+(1:length(agents)), size=n-nrow(insert.dt), replace=TRUE)
    } else {
        agent.order <- rep(1+(1:length(agents)), length.out=n-nrow(insert.dt))
    }
    agent.order <- c(rep(1, nrow(insert.dt)), agent.order)

    list(agents=agents1, order=agent.order)
}
