
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

process_reports <- function(f, nsim, max.iter) {
    done <- integer(length=nsim)
    ks <- vector(mode='list', length=nsim)
    while(!isIncomplete(f)) {
        msg <- readBin(f, 'character')
        msgs <- as.numeric(strsplit(msg, ':', fixed=TRUE)[[1]])
        simi <- msgs[1]
        done[simi] <- msgs[2]
        ks[[simi]] <- c(ks[[simi]], msgs[3])
        if(msgs[2] == max.iter) {
            cat('\r', rep(' ', getOption('width')), sep='', collapse='')
            cat(sprintf('\r(%5.1f%%) %d %s\n',
                        100*sum(done==max.iter) / nsim,
                        simi,
                        agglom(sort(unique(ks[[simi]])))))
            #cat('\r', simi, ' ', agglom(sort(unique(ks[[simi]]))), '\n',
            #    sep='', collapse='')
        }
        curr <- which((done>0) & (done<max.iter))
        cat('\r', paste0(sprintf('%d %6.2f%%', curr, 100*done[curr]/max.iter),
                         collapse=', '),
            sep='')
    }
    cat('\n')
    parallel:::mcexit()
}

make_report_function <- function(f, simi) {
    rf <- function(i, k) {
        writeBin(sprintf('%d:%d:%d', simi, i, k), f)
    }
    rf
}

#' @export
make_categories_runner <- function(positions, nsim=NULL, max.iter=NULL, ...) {
    nsim <- if(is.null(nsim)) max(positions$sim) else nsim
    niter <- if(is.null(max.iter)) max(positions$id) else max.iter

    f <- fifo(tempfile(), open="w+b", blocking=T)
    if (inherits(parallel:::mcfork(), "masterProcess")) {
        process_reports(f, nsim, max.iter=niter)
    }
    result <- parallel::mclapply(1:nsim, function(simi) {
            reportf <- make_report_function(f, simi)
            make_categories(positions=positions[sim==simi],
                            reportf=reportf, ...)
        })
    close(f)
    result
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
make_categories_old <- function(positions,
                                min.iter=positions[agent.id>1, min(id)],
                                max.iter=positions[, max(id)],
                                ic='BIC',
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
        #x <- as.matrix(positions[id <= i, x])
        x <- positions[id <= i, x]

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

#' Make 'optimal' categories for a set of positions
#'
#' Finds the optimal number and position of categories for the evolution of a market
#' Optimal in the sense of information criterion
#' @param positions set of market positions
#' @param init.iter initial iteration to categorize (usually after insertions)
#' @param max.iter maximum iteration to categorize
#' @param package which clustering package to use
#' @param reportf verbosity function (NULL to be quiet)
#' @param ... capture additional arguments
make_categories <- function(positions,
                            min.iter=positions[agent.id>1, min(id)],
                            max.iter=positions[, max(id)],
                            package='mclust',
                            reportf=NULL,
                            ...) {
    # assign NULL mixture to all uncategorized iterations
    #ems <- rep(list(NULL), min.iter-1)
    #ems <- rep(list(NULL), max.iter)
    ems <- vector('list', max.iter)

    for(i in min.iter:max.iter) {
        x <- positions[id <= i, x]
        xcat <- switch(package,
                       'mclust'=categorize.mclust(x, ...),
                       'EMCluster'=categorize.EMCluster(x, ...))

        ems[[i]] <- xcat$mix

        if(!is.null(reportf)) reportf(i, xcat$k)
    }
    ems
}

#' @import mclust
categorize.mclust <- function(x, min.k=1, max.k=9, ...) {
    .min.k <- min.k
    .max.k <- max.k
    res <- mclust::Mclust(x, G=.min.k:.max.k, ...)
    while(res$G == .max.k) {
        .min.k <- .max.k
        .max.k <- 2*.max.k
        res <- mclust::Mclust(x, G=.min.k:.max.k, ...)
    }
    list(mix=res, k=res$G)
}


categorize.EMCluster <- function(x, min.k=1, ic='BIC', ...) {
    nx <- length(x)
    mx <- as.matrix(x)

    catf <- function(k) {
        mix <- EMCluster::init.EM(mx, nclass=k, EMC=EMCluster::.EMC)
        ic <- EMCluster::em.ic(mx, mix)
        list(mix=mix, ic=ic)
    }

    # search space
    k <- min.k
    cata <- catf(k)

    # do while we keep seeing IC improvements
    repeat {
        # make sure we have enough observations to estimate k means and variances (2*k dofs)
        # plus 1 dof for optimization i think? otherwise we have collinearity
        if(2*(k+1) > nx - 1) break
        #if(k+1 > i/2) break

        # search at higher k
        catb <- catf(k+1)

        # get IC differences
        ics <- names(cata$ic)
        ic.diffs <- sapply(ics, function(ic) catb$ic[[ic]] - cata$ic[[ic]])
        names(ic.diffs) <- ics

        # if ic improvement, accept bump in k
        if(ic.diffs[ic] < 0) {
            k <- k+1
            cata <- catb
        } else {
            break
        }
    }
    list(mix=cata$mix, k=k)
}

#' Wrapper around categorization function
#'
#' @param x positions to categorize
#' @param k number of clusters to target
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

get_goms.EMCluster <- function(x, mix, logp=TRUE) {
    res <- sapply(1:mix$nclass, function(ci) {
            cat.ps <- dmixmvn(as.matrix(x), emobj=NULL, log=logp,
                              pi=mix$pi[ci],
                              Mu=mix$Mu[ci,,drop=F],
                              LTSigma=mix$LTSigma[ci,,drop=FALSE])
            #if(peak.difference) {
            #    cat.ps - max(cat.ps)
            #} else {
            #    cat.ps
            #}
            cat.ps
        })
    res
}

get_goms.mclust <- function(x, mix, logp=TRUE) {
    if(logp) {
        log(mix$z)
    } else {
        mix$z
    }
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



