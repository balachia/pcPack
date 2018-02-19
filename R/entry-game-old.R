library(data.table)
library(ggplot2)
library(plyr)
library(purrr)

rm(list=ls())
set.seed(1)

# solver parameters
eps <- 1e-2     # need an epsilon for the root finder
mingap <- 0.3
#mingap <- 0.01
#mingap <- 0

# hyperparameters
# pretty sure b controls risk aversion
a <- 1
b0 <- 1.4
#b0 <- 2
sigma <- 1

c.bs <- list(list(b=1.4, pw=1),
             list(b=2, pw=1))
c.bs <- list(list(b=2, pw=1))
#c.bs <- list(list(b=1.4, pw=1))
bs.string <- sapply(c.bs, function(x) x$b) %>% unique %>% sort %>% paste(collapse='-')

############################################################
# functions

xs <- c(0)
Ws <- c(0)
jumps <- NULL
gaps <- NULL
# find search bounds for root finder
brid.bounds <- function(xl, xr, Wl, Wr, eps0=1e-4, ...) {
    # find bounds for root finder using geometric backoff
    eps <- eps0
    found <- FALSE
    res <- NULL
    while(!found && eps > 1e-10) {
        if(Wr>Wl) {
            rootl <- (xr-xl)/2 + eps
            rootr <- (xr-xl) - eps
        } else {
            rootl <- eps
            rootr <- (xr-xl)/2 - eps
        }
        # check if good, i.e. opposite signed
        critsign <- prod(sign(c(bridcrit(rootl, xl=xl, xh=xr, Wl=Wl, Wh=Wr, ...),
                                bridcrit(rootr, xl=xl, xh=xr, Wl=Wl, Wh=Wr, ...))))
        if(is.na(critsign)) {
            cat(sprintf("missing critsign: [xl %s,xr %s], [l %s, r %s], [Wl %s, Wr %s] ->\n\t(%s, %s)\n", xl, xr, rootl, rootr, Wl, Wr,
                        bridcrit(rootl, xl=xl, xh=xr, Wl=Wl, Wh=Wr, ...),
                        bridcrit(rootr, xl=xl, xh=xr, Wl=Wl, Wh=Wr, ...)))
            stop("Can\'t evaluate bridge criterion")
        }
        #if(prod(sign(c(bridcrit(rootl, xl=xl, xh=xr, Wl=Wl, Wh=Wr, ...), bridcrit(rootr, xl=xl, xh=xr, Wl=Wl, Wh=Wr, ...)))) < 0) {
        if(critsign < 0) {
            res <- list()
            res$l <- rootl
            res$r <- rootr
            found <- TRUE
        } else {
            #cat(sprintf("reducing eps: %s -> %s\n", eps, eps/2))
            eps <- eps/2
        }
    }
    res
}

# process information at new split
# i - curr index
# splitidx - index of best jump
# Wnew - value at new jump
split.interval <- function(i, splitidx, xnew, Wnew, bnew, info) {
    # split old interval into left part (splitidx) and right part (i)
    info.new <- data.table(id=c(i, splitidx),
                           xl=c(xnew, info[splitidx, xl]),
                           xr=c(info[splitidx, xr], xnew),
                           Wl=c(Wnew, info[splitidx, Wl]),
                           Wr=c(info[splitidx, Wr], Wnew),
                           bl=c(bnew, info[splitidx, bl]),
                           br=c(info[splitidx, br], bnew),
                           key="id")
    info.new
}

# reconstruct expectations at new split
split.expect <- function(i, splitidx, info, bs, open.min=5e-4, open.max=1e5) {
    plans.new <- lapply(bs, function(bx) {
        b <- bx$b
        plan.new <- data.table(id=c(i, splitidx), key="id")

        # get new interval expectations for the lefter interval (splitidx)
        if (is.infinite(info[J(splitidx), xl])) {
            tryCatch({
                d0 <- uniroot(opencrit, c(open.min, open.max), W0=info[J(splitidx), Wr], b=b)$root
            }, error=function(e) {
                cat(sprintf('bad crit sign (W0: %s): %s, %s',
                            info[J(splitidx), Wr],
                            opencrit(open.min, W0=info[J(splitidx), Wr], b=b),
                            opencrit(open.max, W0=info[J(splitidx), Wr], b=b)))
                stop(e)
            })
            #d0 <- uniroot(opencrit, c(5e-3, 1e2), W0=info[J(splitidx), Wr], b=b)$root
            Eu0 <- Euopen(d0, info[J(splitidx), Wr], b=b)
            #info.new[J(splitidx), `:=`(delta=d0, Eu=Eu0, xd=xr-d0)]
            plan.new[J(splitidx), `:=`(delta=d0, Eu=Eu0, xd=info[J(splitidx), xr]-d0)]
        } else {
            res <- tryCatch({
                root <- brid.bounds(xl=info[J(splitidx), xl], xr=info[J(splitidx), xr],
                                    Wl=info[J(splitidx), Wl], Wr=info[J(splitidx), Wr], b=b)
                gaps <- c(gaps, (info[J(splitidx), xr-xl]))
                d0 <- uniroot(bridcrit, c(root$l, root$r), b=b,
                              xl=info[J(splitidx), xl], xh=info[J(splitidx), xr],
                              Wl=info[J(splitidx), Wl], Wh=info[J(splitidx), Wr])$root
                Eu0 <- Eubrid(d0, info[J(splitidx), xl], info[J(splitidx), xr],
                              info[J(splitidx), Wl], info[J(splitidx), Wr], b=b)
                list(d0, Eu0)
            }, error=function(e) {
                d0 <- info[J(splitidx), (xr-xl)/2]
                Eu0 <- -Inf
                list(d0, Eu0)
            }, finally={})
            d0 <- res[[1]]
            Eu0 <- res[[2]]
            #if(info[J(splitidx), xr-xl>mingap]){
            #    root <- brid.bounds(xl=info[J(splitidx), xl], xr=info[J(splitidx), xr],
            #                        Wl=info[J(splitidx), Wl], Wr=info[J(splitidx), Wr], b=b)
            #    gaps <- c(gaps, (info[J(splitidx), xr-xl]))
            #    d0 <- uniroot(bridcrit, c(root$l, root$r), b=b,
            #                  xl=info[J(splitidx), xl], xh=info[J(splitidx), xr],
            #                  Wl=info[J(splitidx), Wl], Wh=info[J(splitidx), Wr])$root
            #    Eu0 <- Eubrid(d0, info[J(splitidx), xl], info[J(splitidx), xr],
            #                  info[J(splitidx), Wl], info[J(splitidx), Wr], b=b)
            #} else {
            #    d0 <- info[J(splitidx), (xr-xl)/2]
            #    Eu0 <- -Inf
            #}
            plan.new[J(splitidx), `:=`(delta=d0, Eu=Eu0, xd=info[J(splitidx), xl]+d0)]
        }

        # get new interval expectations for the righter interval (i)
        if (is.infinite(info[J(i), xr])) {
            tryCatch({
                d0 <- uniroot(opencrit, c(open.min, open.max), W0=info[J(i), Wl], b=b)$root
            }, error=function(e) {
                cat(sprintf('bad crit sign (W0: %s): %s, %s',
                            info[J(i), Wl],
                            opencrit(open.min, W0=info[J(i), Wl], b=b),
                            opencrit(open.max, W0=info[J(i), Wl], b=b)))
                stop(e)
            })
            #d0 <- uniroot(opencrit, c(5e-2, 1e2), W0=info[J(i), Wl], b=b)$root
            Eu0 <- Euopen(d0, info[J(i), Wl], b=b)
            plan.new[J(i), `:=`(delta=d0, Eu=Eu0, xd=info[J(i), xl]+d0)]
        } else {
            res <- tryCatch({
                root <- brid.bounds(xl=info[J(i), xl], xr=info[J(i), xr],
                                    Wl=info[J(i), Wl], Wr=info[J(i), Wr], b=b)
                gaps <- c(gaps, (info[J(i), xr-xl]))
                d0 <- uniroot(bridcrit, c(root$l, root$r), b=b,
                              xl=info[J(i), xl], xh=info[J(i), xr],
                              Wl=info[J(i), Wl], Wh=info[J(i), Wr])$root
                Eu0 <- Eubrid(d0, info[J(i), xl], info[J(i), xr],
                              info[J(i), Wl], info[J(i), Wr], b=b)
                list(d0, Eu0)
            }, error=function(e) {
                d0 <- info[J(i), (xr-xl)/2]
                Eu0 <- -Inf
                list(d0, Eu0)
            }, finally={})
            d0 <- res[[1]]
            Eu0 <- res[[2]]
            #if(info[J(i), xr-xl>mingap]){
            #    root <- brid.bounds(xl=info[J(i), xl], xr=info[J(i), xr],
            #                        Wl=info[J(i), Wl], Wr=info[J(i), Wr], b=b)
            #    gaps <- c(gaps, (info[J(i), xr-xl]))
            #    d0 <- uniroot(bridcrit, c(root$l, root$r), b=b,
            #                  xl=info[J(i), xl], xh=info[J(i), xr],
            #                  Wl=info[J(i), Wl], Wh=info[J(i), Wr])$root
            #    Eu0 <- Eubrid(d0, info[J(i), xl], info[J(i), xr],
            #                  info[J(i), Wl], info[J(i), Wr], b=b)
            #} else {
            #    d0 <- info[J(i), (xr-xl)/2]
            #    Eu0 <- -Inf
            #}
            plan.new[J(i), `:=`(delta=d0, Eu=Eu0, xd=info[J(i), xl]+d0)]
        }
        plan.new
    })
    plans.new
}

make.split <- function(newidx, splitidx, x, W, b, info, bs) {
    info.new <- split.interval(newidx, splitidx, x, W, b, info)
    info[info.new, `:=`(xl=i.xl, xr=i.xr, Wl=i.Wl, Wr=i.Wr, bl=i.bl, br=i.br)]

    # recalculate agent expectations
    plans.new <- split.expect(newidx, splitidx, info, bs)
    for(bi in 1:length(bs)) {
        bs[[bi]]$plan[plans.new[[bi]], `:=`(delta=i.delta, Eu=i.Eu, xd=i.xd)]
    }

    list(info=info, bs=bs)
}

# find insertion point
insert <- function(x, W, info) {
    res <- list(idx=info[(is.finite(xr) | is.finite(xl)) & x > xl & x < xr, id], x=x, W=W)
    if(length(res$idx)==0) res$idx <- 1
    res
}

best <- function(plan, info) {
    #bestidx <- info[, order(Eu, decreasing=TRUE)[1]]
    bestidx <- plan[, order(Eu, decreasing=TRUE)[1]]
    if(is.infinite(info[bestidx, xl])) {
        Wnew <- info[bestidx, Wr] + openjump(plan[bestidx, delta])
    } else if (is.infinite(info[bestidx, xr])) {
        Wnew <- info[bestidx, Wl] + openjump(plan[bestidx, delta])
    } else {
        Wnew <- bridjump(plan[bestidx, delta],
                         info[bestidx, xl], info[bestidx, xr],
                         info[bestidx, Wl], info[bestidx, Wr])
    }
    list(idx=bestidx, x=plan[bestidx, xd], W=Wnew)
}

# insertions
init.gap <- 4
insert.xs <- c(0, 50, -init.gap, init.gap, 50-init.gap, 50+init.gap)
insert.Ws <- c(0, 0, rep(-sqrt(init.gap)/2, 4))
insert.dt <- data.table(id=1:length(insert.xs) + 1,
                        x=insert.xs,
                        W=insert.Ws,
                        b=as.numeric(NA))

simrun <- function(n, bs=c.bs) {
    # initialize information table
    info <- data.table(id=1:(n+1), xl=-Inf, xr=Inf,
                       Wl=as.numeric(NA), Wr=as.numeric(NA),
                       bl=as.numeric(NA), br=as.numeric(NA))

    setkey(info, id)
    d0 <- uniroot(opencrit, c(5e-2, 1e2), W0=0)$root
    Eu0 <- Euopen(d0, 0)
    #info[1, `:=`(xr=0, Wr=0, delta=d0, Eu=Eu0, xd=-d0)]
    #info[2, `:=`(xl=0, Wl=0, delta=d0, Eu=Eu0, xd=d0)]
    #info[1, `:=`(xr=0, Wr=0)]
    #info[2, `:=`(xl=0, Wl=0)]

    # create expectation for all agent types
    #bs <- list(list(b=1.4, pw=1),
    #           list(b=2, pw=1))
    bprobs <- sapply(bs, function (x) x$pw)

    # augment bs with expectation table
    bs <- lapply(bs, function (bx) {
            b <- bx$b
            d0 <- uniroot(opencrit, c(5e-2, 1e2), W0=0, b=b)$root
            Eu0 <- Euopen(d0, 0, b=b)
            plan <- data.table(id=1:(n+2))
            setkey(plan, id)
            plan[1, `:=`(delta=d0, Eu=Eu0, xd=-d0)]
            plan[2, `:=`(delta=d0, Eu=Eu0, xd=d0)]
            bx$plan <- plan
            bx
        })

    agg.dists <- NULL

    # simulate
    #cat('\n')
    for (i in (1:n + 1)) {
        #plot(c(info[is.finite(xr), xr], info[!is.na(xd), xd]),
             #c(info[is.finite(xr), Wr], info[!is.na(xd), Eu]),
             #pch=c(rep(15, i-2), rep(20, i-1)), ylim=c(-45, 50))
        #plot(info[is.finite(xr), xr],
             #info[is.finite(xr), Wr],
             #pch=rep(15, i-2), ylim=c(-45, 50))

        # pick agent
        bidx <- sample.int(length(bs), size=1, prob=bprobs)
        b <- bs[[bidx]]$b
        plan <- bs[[bidx]]$plan

        cat(sprintf("\r[%s] %6.2f %% (bidx %d -- %3.1f)\r[%s\r",
                    paste0(rep("-", 50), collapse=""),
                    (i-1)/n*100,
                    bidx, b,
                    paste0(rep("=", floor(50*(i-1)/n)), collapse="")))

        # pick best option, and split
        #bestidx <- info[, order(Eu, decreasing=TRUE)[1]]
        #bestidx <- plan[, order(Eu, decreasing=TRUE)[1]]
        #if(is.infinite(info[bestidx, xl])) {
            ##Wnew <- info[bestidx, Wr] + openjump(info[bestidx, delta])
            #Wnew <- info[bestidx, Wr] + openjump(plan[bestidx, delta])
        #} else if (is.infinite(info[bestidx, xr])) {
            ##Wnew <- info[bestidx, Wl] + openjump(info[bestidx, delta])
            #Wnew <- info[bestidx, Wl] + openjump(plan[bestidx, delta])
        #} else {
            ##Wnew <- bridjump(info[bestidx, delta],
                             ##info[bestidx, xl], info[bestidx, xr],
                             ##info[bestidx, Wl], info[bestidx, Wr])
            #Wnew <- bridjump(plan[bestidx, delta],
                             #info[bestidx, xl], info[bestidx, xr],
                             #info[bestidx, Wl], info[bestidx, Wr])
        #}

        # are we inserting at this iteration?
        if(nrow(insert.dt[id==i])>0) {
            split.at <- insert(insert.dt[id==i, x], insert.dt[id==i, W], info)
            split.at$b <- insert.dt[id==i, b]
        } else {
            split.at <- best(plan, info)
            split.at$b <- b
        }

        #res <- make.split(i, bests$idx, bests$x, bests$W, b, info, bs)
        res <- make.split(i, split.at$idx, split.at$x, split.at$W, split.at$b, info, bs)
        info <- res$info
        bs <- res$bs

        # make dataset centered on points, and aggregate it over runs
        dtl <- info[is.finite(xl)]
        dtr <- info[is.finite(xr)]
        setnames(dtl, c("xl", "Wl", "bl"), c("x", "W", "b"))
        setnames(dtr, c("xr", "Wr", "br"), c("x", "W", "b"))
        setkey(dtl, x)
        setkey(dtr, x)
        dist.dt <- dtl[dtr, list(iter=i-1, x, xl, xr, W, b, dist=pmin(xr-x, x-xl))]
        dist.dt[, ct := ct(dist)]
        agg.dists <- rbind(agg.dists, dist.dt)
        #agg.dists <- rbind(agg.dists, info[is.finite(xr-xl), list(iter=i-2, xl, Wl, bl, dist=xr-xl, ct=ct(xr-xl))])

        ## make split
        #info.new <- split.interval(i, bestidx, plan[bestidx, xd], Wnew, info)
        #info[info.new, `:=`(xl=i.xl, xr=i.xr, Wl=i.Wl, Wr=i.Wr)]

        ## recalculate agent expectations
        #plans.new <- split.expect(i, bestidx, info, bs)
        #for(bi in 1:length(bs)) {
            #bs[[bi]]$plan[plans.new[[bi]], `:=`(delta=i.delta, Eu=i.Eu, xd=i.xd)]
        #}
    }
    cat("\n")

    agg.dists[, V:=W+ct]
    list(info=info, agg.dists=agg.dists)
}

############################################################
# script

niter <- 1e2
niter <- 1e3
nsim <- 1e2

#ress <- lapply(1:nsim, function(simiter) {
#        cat(sprintf("%d/%d\n", simiter, nsim))
#        res <- simrun(niter)
#        res$info[, simiter := simiter]
#        res$agg.dists[, simiter := simiter]
#        res
#    })

#res <- list(info=rbindlist(lapply(ress, function (x) x$info)),
#            agg.dists=rbindlist(lapply(ress, function(x) x$agg.dists)))

## make category like things
#cat.thresh <- 1.5
##res$agg.dists[order(x), catgroup:=cumsum((x-xl)>cat.thresh), by=list(simiter, iter)]
#res$agg.dists[order(simiter, iter, x), catgroup:=cumsum( (x-xl)>cat.thresh)]
#res$agg.dists[, `:=`(catpeak.x=x[which.max(W)],
#                    catpeak.W=W[which.max(W)])
#              , by=catgroup]

##saveRDS(res, sprintf("Rds/brownian-results-%dsim-%diter.Rds", nsim, niter))
#saveRDS(res, sprintf("Rds/brownian-results-%dsim-%diter-b%s.Rds", nsim, niter, bs.string))
