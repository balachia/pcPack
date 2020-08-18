############################################################
# make jsonable

#' @export
jsonable <- function(...) {
    UseMethod('jsonable')
}

#' @export
jsonable.default <- function(x, ...) x

#' @export
jsonable.simulation_listing <- function(x, ...) {
    x$.class <- class(x)
    class(x) <- 'list'
    x$packageVersion <- as.character(x$packageVersion)
    x$agents <- lapply(x$agents, jsonable)
    x$categorizations <- lapply(x$categorizations, jsonable)
    x
}

#' @export
jsonable.categorization_listing <- function(x, ...) {
    x$.class <- class(x)
    class(x) <- 'list'
    x$packageVersion <- as.character(x$packageVersion)
    x
}

#' @export
jsonable.agent <- function(x, ...) {
    x$.class <- class(x)
    class(x) <- 'list'
    x
}

############################################################
# dejsonize

dejsonize <- function(x, ...) {
    
}


############################################################
# one-time converter

parameters2json <- function(parameters_file=default_parameters_file, ...) {
    params <- read_parameters(file=parameters_file)
    params <- lapply(params, jsonable)
    print(params)

    write(jsonlite::toJSON(params, pretty=TRUE), 'parameters.json')

}
