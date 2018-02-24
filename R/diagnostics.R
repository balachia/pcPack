
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
