
#' pairwise distance matrix for vectors
vectorized_pdist <- function(A,B) {
    an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
    bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
 
    m = nrow(A)
    n = nrow(B)
 
    tmp = matrix(rep(an, n), nrow=m) 
    tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
    sqrt( tmp - 2 * tcrossprod(A,B) )
}

#' pairwise distance matrix for scalars
pdist <- function(x,y) {
    m = length(x)
    n = length(y)

    tmp = matrix(rep(x^2, n), nrow=m) 
    tmp = tmp +  matrix(rep(y^2, m), nrow=m, byrow=TRUE)
    sqrt( tmp - 2 * tcrossprod(x,y) )
}

#' distance to closest neighbor
#'
#' @export
#' @importFrom matrixStats rowOrderStats
neighbor_distance <- function(x, min.n=1, max.n=Inf) {
    ds <- pdist(x,x)
    diag(ds) <- Inf
    lapply(min.n:min(max.n, length(x)), function(n) {
        rowOrderStats(ds, 1:n, 1:n, which=1)
    })
}

#' gom category statistics
#'
#' @export
gom_statistics <- function(goms) {
    stats <- data.table(gom1=goms[,1],
                        gom2=goms[,2])
    stats[, `:=`(miscat=exp(gom2) - exp(gom1),
                 spanner=exp(gom2-gom1))]
}

#' peak category grade of membership
#'
#' @export
peak_gom <- function(x) {
    x[,1]
}

#' second category grade of membership
#'
#' @export
second_gom <- function(x) {
    x[,2]
}


#' spanner function
#'
#' @export
spanner <- function(x) {
    x[,1] - x[,2]
}

#' expand finals positions to all iterations
#'
#' @export
expand_positions <- function(positions, sims, iters) {
    expander <- CJ(sim=1:sims, iter=rep(1:iters, 1:iters))
    expander[, id := 1:.N, by=list(sim, iter)]

    big.positions <- expander[positions, on=c('sim', 'id')]
    setkey(big.positions, 'sim', 'iter')
    big.positions
}

#' add category definitions to positions table
#'
#' @export
position_statistics <- function(positions, goms, min.iter = NULL) {
    if(is.null(min.iter)) {
        min.iter <- sum(sapply(goms[[1]], is.null)) + 1
    }
    if(min.iter > 1) { goms <- lapply(goms, function(sim.goms) tail(sim.goms, -min.iter + 1)) }

    big.goms <- do.call(rbind,
                        lapply(goms, function(sim.goms) do.call(rbind,
                                                                lapply(sim.goms, function(x) x$goms))))
    big.ords <- do.call(rbind,
                        lapply(goms, function(sim.ords) do.call(rbind,
                                                                lapply(sim.ords, function(x) x$ords))))

    gom_stats <- gom_statistics(big.goms)
    positions[iter >= min.iter, names(gom_stats) := gom_stats]
    positions[iter >= min.iter, cat.id := big.ords[,1]]
    #positions[, dW := W - mean(W), by=list(sim, iter, cat.id)]

    # binarized
    positions[, b.miscat := miscat > -0.5]
    positions[, b.spanner := spanner > 0.01]

    # split out to design matrix
    positions[, b.pure := (!b.miscat)]
    positions[, b.span := (b.miscat & b.spanner)]
    positions[, b.none := (b.miscat & !b.spanner)]

    # factor
    positions[ (b.pure), cat.type := 'pure']
    positions[ (b.none), cat.type := 'none']
    positions[ (b.span), cat.type := 'span']

    positions
}

