default_data_dir <- '../data/'
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

    cat(sprintf('\tcats:'))
    if(length(x$categorizations) > 0) cat('\r') else cat('\n')
    for(catl in x$categorizations) { cat('\t\t'); print(catl) }

    #cat(sprintf('\t# cats:\t%d\n', length(x$categorizations)))
}

#' @export
make_categorization_listing <- function(simulation_digest, package,
                                        overwrite=FALSE,
                                        ic=NA,
                                        min.iter=NA, max.iter=NA,
                                        min.k=1, max.k=Inf,
                                        ...) {
    res <- list()
    class(res) <- 'categorization_listing'

    res$packageVersion <- packageVersion('prodcatpack')
    res$simulation_digest <- simulation_digest
    res$package <- package
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
print.categorization_listing <- function(x, full=FALSE, ...) {
    if(full) cat(sprintf('%s (%s) ', x$digest, x$simulation_digest))
    cat(sprintf('%s', x$package))
    if(!is.na(x$ic)) cat(sprintf(', ic %s', x$ic))
    if(!is.na(x$min.iter) || !is.na(x$max.iter)) {
        cat(', iter ')
        if(!is.na(x$min.iter)) cat(sprintf('%.0f.', x$min.iter))
        cat('.')
        if(!is.na(x$max.iter)) cat(sprintf('.%.0f', x$max.iter))
    }
    if((x$min.k > 1) || is.finite(x$max.k)) {
        cat(', k ')
        if(x$min.k > 1) cat(sprintf('%.0f.', x$min.k))
        cat('.')
        if(is.finite(x$max.k)) cat(sprintf('.%.0f', x$max.k))
    }
    cat('\n')
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
                                    overwrite=FALSE,
                                    ...) {
    parameters <- read_parameters(parameters_file)
    cat_listing <- make_categorization_listing(simulation_digest=simulation_digest, ...)
    if( overwrite || is.null(parameters[[simulation_digest]]$categorizations[[cat_listing$digest]]) ) {
        parameters[[simulation_digest]]$categorizations[[cat_listing$digest]] <- cat_listing
    }
    #parameters[[simulation_digest]]$categorizations <- c(parameters[[simulation_digest]]$categorizations, list(cat_listing))

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

category_to_simulation_map <- function(parameters_file=default_parameters_file) {
    params <- read_parameters(file=parameters_file)
    cats <- list()
    for(siml in params) {
        for(catl in siml$categorizations) cats[[catl$digest]] <- siml$digest
    }
    cats
}

#' @export
load_simulation_data <- function(digest,
                                 parameters_file=default_parameters_file,
                                 data_dir=default_data_dir) {
    siml <- read_parameters(file=parameters_file)[[digest]]
    dat <- readRDS(sprintf('%s/simulations/%s.Rds', data_dir, digest))
    list(listing=siml, data=dat)
}

#' @export
load_categorization_data <- function(digest,
                                 parameters_file=default_parameters_file,
                                 data_dir=default_data_dir) {
    sim_digest <- category_to_simulation_map(parameters_file)[[digest]]
    catl <- read_parameters(file=parameters_file)[[sim_digest]]$categorizations[[digest]]
    dat <- readRDS(sprintf('%s/categorizations/%s.Rds', data_dir, digest))
    sim_dat <- load_simulation_data(sim_digest, parameters_file=parameters_file, data_dir=data_dir)
    list(listing=catl, data=dat, simulation=sim_dat)
}
