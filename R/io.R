default_parameters_file <- '../data/parameters.Rds'

#' @export
read_parameters <- function(file=default_parameters_file) {
    readRDS(file)
}

#' @export
write_parameters <- function(parameters, file=default_parameters_file) {
    saveRDS(parameters, file)
}

make_parameters <- function(file=default_parameters_file) {
    if(!dir.exists(dirname(file))) dir.create(dirname(file))
    write_parameters(list(), file)
}


#' @export
make_simulation_listing <- function(nsim, niter, peaks, seed, agents, ...) {
    res <- list()
    class(res) <- 'simulation_listing'

    res$packageVersion <- packageVersion('prodcatpack')
    res$nsim <- nsim
    res$niter <- niter
    res$peaks <- peaks
    res$seed <- seed
    res$agents <- agents

    # make digest
    res$digest <- digest::digest(res)

    # attach empty categorization registry
    res$categorizations <- list()

    res
}

#' @export
print.simulation_listing <- function(x, ...) {
    cat(sprintf('%s\n', x$digest))
    cat(sprintf('\tversion:\t%s\n', x$packageVersion))
    cat(sprintf('\tsims:\t%d\n', x$nsim))
    cat(sprintf('\titers:\t%d\n', x$niter))
    cat(sprintf('\tpeaks:\t%d\n', x$peaks))
    cat(sprintf('\tseed:\t%d\n', x$seed))
    cat(sprintf('\tagents:\t'))
    if(!is.null(x$agents)) {
        cat('\r')
        for(ag in x$agents) { cat('\t\t'); print(ag) }
    } else {
        cat(sprintf('NA\n'))
    }

    cat(sprintf('\t# cats:\t%d\n', length(x$categorizations)))
}

#' @export
make_categorization_listing <- function(simulation_digest, ic,
                                        min.iter=NA, max.iter=NA,
                                        min.k=1, max.k=Inf,
                                        ...) {
    res <- list()
    class(res) <- 'categorization_listing'

    res$packageVersion <- packageVersion('prodcatpack')
    res$simulation_digest <- simulation_digest
    res$ic <- ic
    res$min.iter <- min.iter
    res$max.iter <- max.iter
    res$min.k <- min.k
    res$max.k <- max.k

    # make digest
    res$digest <- digest::digest(res)

    res
}

#' @export
register_simulation <- function(parameters_file=default_parameters_file,
                                overwrite=FALSE,
                                ...) {
    parameters <- read_parameters(parameters_file)
    sim_listing <- make_simulation_listing(...)

    # if sim doesn't already exist
    if(is.null(parameters[[sim_listing$digest]])) {
        #parameters <- c(parameters, list(sim_listing))
        #names(parameters)[length(parameters)] <- sim_listing$digest
        parameters[[sim_listing$digest]] <- sim_listing
    } else {
        if(overwrite) {
            parameters[[sim_listing$digest]] <- sim_listing
        } else {
            # else ignore new simulation
            warning(sprintf('Simulation (%s) already exists in parameters file.',
                            sim_listing$digest))
        }
    }

    write_parameters(parameters, file=parameters_file)
    sim_listing
}

#' @export
register_categorization <- function(simulation_digest,
                                    parameters_file=default_parameters_file,
                                    ...) {
    parameters <- read_parameters(parameters_file)
    cat_listing <- make_categorization_listing(simulation_digest=simulation_digest, ...)
    parameters[[simulation_digest]]$categorizations <- c(parameters[[simulation_digest]]$categorizations, list(cat_listing))

    write_parameters(parameters, file=parameters_file)
    cat_listing
}

#' @export
#find_simulation <- function(parameters_file=default_parameters_file, ...) {
find_simulation <- function(parameters, ...) {
    #parameters <- read_parameters(parameters_file)
    l <- list(...)
    
    # search parameters
    for(i in seq_along(l)) {
        pname <- names(l)[i]
        pval <- l[[i]]
        parameters <- Filter(function(x) { x[[pname]] == pval }, parameters)
    }

    parameters
}

#' @export
find_simulation_by_standard_agent <- function(parameters, agents) {
    `%subset%` <- function(x, y) all(x %in% y)
    `%set.equal%` <- function(x, y) x %subset% y && y %subset% x

    parameters <- Filter(function(p) p$agents %set.equal% agents, parameters)
    parameters
}

#' @export
deregister_simulations <- function(digests, parameters_file=default_parameters_file, ...) {
    parameters <- read_parameters(parameters_file)

    keep.names <- Filter(function(name) !(name %in% digests), names(parameters))
    keep.parameters <- parameters[keep.names]

    write_parameters(keep.parameters, file=parameters_file)
}
