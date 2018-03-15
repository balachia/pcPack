
agglom.reduc <- function(xs,x) {
    idx <- length(xs)
    if(idx == 0) {
        xs <- list(c(x,x))
    } else {
        if((x - xs[[idx]][2]) == 1) {
            xs[[idx]] = c(xs[[idx]][1],x)
        } else {
            xs[[idx+1]] = c(x,x)
        }
    }
    xs
}

agglom <- function(xs) {
    res <- Reduce(agglom.reduc, xs, list())
    res <- lapply(res, function(x) { if(x[1]==x[2]) { x[1] } else { paste0(x[1],'-',x[2]) }})
    paste0(res, collapse=',')
}

#' Make 'optimal' categories for a set of positions
#'
#' Finds the optimal number and position of categories for the evolution of a market
#' Optimal in the sense of information criterion
#' @param positions set of market positions
#' @param init.iter initial iteration to categorize (usually after insertions)
#' @param max.iter maximum iteration to categorize
#' @param method information criterion method
#' @param ... captire additional arguments
#' @export
make_categories <- function(positions,
                            min.iter=positions[agent.id>1, min(id)],
                            max.iter=positions[, max(id)],
                            ic='AIC',
                            verbose.prefix='',
                            ...) {
    # verbosity calculations
    base.line.format <- '\r%sITER %%%ds :: k %%-%ds'
    max.pad <- 2
    line.format <- sprintf(base.line.format, verbose.prefix, 1+floor(log10(max.iter)), max.pad)

    # assign NULL mixture to all uncategorized iterations
    ems <- rep(list(NULL), min.iter-1)
    ks <- NULL
    #ems <- rep(list(NULL), max.iter)

    for(i in min.iter:max.iter) {
        x <- as.matrix(positions[id <= i, x])

        #cat('\n')

        # search space
        k <- 1
        cat(sprintf(line.format, i, agglom(c(ks, k))))
        #mixa <- EMCluster::init.EM(x, nclass=k, EMC=EMCluster::.EMC)
        #ic.a <- EMCluster::em.ic(x,mixa)
        cata <- categorize.em(x, k)

        # do while we keep seeing IC improvements
        repeat {
            if(k+1 > i/2) break

            cat(sprintf(line.format, i, agglom(c(ks, k))))
            # search at higher k
            #mixb <- EMCluster::init.EM(x, nclass=k+1, EMC=EMCluster::.EMC)
            #ic.b <- EMCluster::em.ic(x,mixb)
            catb <- categorize.em(x, k+1)

            # get IC differences
            ics <- names(cata$ic)
            ic.diffs <- sapply(ics, function(ic) catb$ic[[ic]] - cata$ic[[ic]])
            names(ic.diffs) <- ics

            if(ic.diffs[ic] < 0) {
                k <- k+1
                cata <- catb
                #mixa <- mixb
                #ic.a <- ic.b
            } else {
                break
            }
        }

        ems <- c(ems, list(cata$mix))
        ks <- sort(unique(c(ks, k)))
        if(nchar(agglom(ks)) + 3 > max.pad) max.pad <- nchar(agglom(ks)) + 3
        line.format <- sprintf(base.line.format, verbose.prefix, 1+floor(log10(max.iter)), max.pad)
        #ems[[i]] <- list(cata$mixa)
    }
    ems
}

#' Wrapper around categorization function
#'
#' @param x positions to categorize
#' @param k number of clusters to target
#' @importFrom EMCluster init.EM em.ic
categorize.em <- function(x, k) {
    mix <- EMCluster::init.EM(x, nclass=k, EMC=EMCluster::.EMC)
    ic <- EMCluster::em.ic(x, mix)
    list(mix=mix, ic=ic)
}

#' Assigns category grade of membership to positions
#'
#' @param x positions
#' @param mix category mixture distribution
#' @param peak.difference grade of membership = difference from peak GOM within category
#' @importFrom EMCluster dmixmvn
#' @export
categorize_positions <- function(x, mix, peak.difference=TRUE) {
    res <- sapply(1:mix$nclass, function(ci) {
            cat.ps <- dmixmvn(as.matrix(x), emobj=NULL, log=TRUE,
                              pi=mix$pi[ci],
                              Mu=mix$Mu[ci,,drop=F],
                              LTSigma=mix$LTSigma[ci,,drop=FALSE])
            if(peak.difference) {
                cat.ps - max(cat.ps)
            } else {
                cat.ps
            }
        })
    res
}

#' Find GOM for two two categories
#'
#' @param goms grades of membership for all categories
#' @return list(goms=GOM for top, second categories; ords=category ids)
#' @export
top_two_categories <- function(goms) {
    # find distance to second category
    # within each row, sort goms from highest to lowest
    goms_ords <- sapply(1:nrow(goms), function(ri) {
            ord <- order(goms[ri,], decreasing=TRUE)
            c(goms[ri,ord], ord) }
        )
    gom_ords <- t(gom_ords)

    #if(iter.cats$nclass > 1) {
    if(ncol(goms) > 1) {
        #res <- goms_ords[,c(1:2, iter.cats$nclass + 1:2)]
        #res <- goms_ords[,c(1:2, ncol(goms) + 1:2)]
        res <- list(goms=goms_ords[,1:2], ords=goms_ords[,1:2 + ncol(goms)])
    } else {
        # if we only have one category, second category has -Inf GOM, and no order
        #res <- cbind(goms[,1], -Inf, 1, NA)
        res <- list(goms=cbind(goms[,1], -Inf), ords=cbind(goms_ords[,2], NA))
    }
    res
}



