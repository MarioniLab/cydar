# This tests the normalization machinery in normalizeBatch.
# library(cydar); library(testthat); source("test-norm.R")

# Mocking up some data.
nmarkers <- 10
marker.names <- paste0("X", seq_len(nmarkers))
all.x <- list()

for (b in paste0("Batch", 1:3)) { # 3 batches
    nsamples <- 10
    sample.names <- paste0("Y", seq_len(nsamples))
    trans.shift <- runif(nmarkers, 0, 1)
    trans.grad <- runif(nmarkers, 1, 2)
    x <- list()
    for (i in sample.names) {
        ncells <- round(runif(1, 500, 2000))
        ex <- matrix(rgamma(nmarkers*ncells, 2, 2), nrow=nmarkers)
        ex <- t(ex*trans.grad + trans.shift)
        colnames(ex) <- marker.names
        x[[i]] <- ex
    }   
    all.x[[b]] <- x
}

# Each batch contains different composition/ordering of groups
batch.comp <- list( 
    factor(rep(1:2, c(3,7))),
    factor(rep(1:2, c(7,3))),
    factor(rep(1:2, 5))
)

# Common values to be used in downstream tests.
batch.out <- lapply(all.x, cydar:::.pull_out_data)
test.wts <- cydar:::.computeCellWeights(batch.out, batch.comp)

###########################################################

test_that("Normalization weights are correctly calculated", {
    getTargetComposition <- function(batch.comp) {
        tab <- table(rep(seq_along(batch.comp), lengths(batch.comp)), unlist(batch.comp))
        tab[,colSums(tab > 0)!=nrow(tab)] <- 0
        colMeans(tab)
    }

    averages <- getTargetComposition(batch.comp)
    for (b in seq_along(batch.comp)) {
        ncells <- unname(sapply(all.x[[b]], nrow))
        cur.b <- batch.comp[[b]]
        cur.n <- unname(averages[cur.b]/as.vector(table(cur.b)[cur.b]))
        expect_equal(cur.n/ncells, test.wts[[b]])
    }

    # Now with a missing group in one batch.
    missing.batch.comp <- list( 
        factor(rep(1:2, c(2,8))),
        factor(rep(1:2, c(6,4))),
        factor(rep(1, 10), levels=1:2)
    )
    test.wts <- cydar:::.computeCellWeights(batch.out, missing.batch.comp)

    averages <- getTargetComposition(missing.batch.comp)
    for (b in seq_along(missing.batch.comp)) { 
        ncells <- unname(sapply(all.x[[b]], nrow))
        cur.b <- missing.batch.comp[[b]]
        cur.n <- unname(averages[cur.b]/as.vector(table(cur.b)[cur.b]))
        cur.n[cur.n==0] <- NA
        expect_equal(cur.n/ncells, test.wts[[b]])
    }

    # Now with a unique group in one batch.
    unique.batch.comp <- list( 
        factor(rep(1:2, c(6, 4)), levels=1:3),
        factor(rep(1:2, c(9, 1)), levels=1:3),
        factor(rep(c(1,3), c(6, 4)), levels=1:3)
    )
    test.wts <- cydar:::.computeCellWeights(batch.out, unique.batch.comp)

    averages <- getTargetComposition(unique.batch.comp)
    for (b in 1:3) {
        ncells <- unname(sapply(all.x[[b]], nrow))
        cur.b <- unique.batch.comp[[b]]
        cur.n <- unname(averages[cur.b]/as.vector(table(cur.b)[cur.b]))
        cur.n[cur.n==0] <- NA
        expect_equal(cur.n/ncells, test.wts[[b]])
    } 
})


test_that("Marker data extraction works as expected", {
    extractMarker <- function(exprs, m) {
        unlist(lapply(exprs, FUN=function(x) { x[,m] }), use.names=FALSE)
    }

    batch.weights <- list(runif(10), runif(10), runif(10))
    FUN <- cydar:::.pullOutMarkers(batch.out, batch.weights)

    ref.weights <- batch.weights
    for (i in seq_along(ref.weights)) {
        ncells <- vapply(batch.out[[i]]$exprs, nrow, 0L)
        ref.weights[[i]] <- rep(batch.weights[[i]], ncells)
    }

    for (m in marker.names) {
        test <- FUN(m)
        expect_identical(ref.weights, test$weights)
        for (i in seq_along(batch.out)) {
            expect_identical(test$exprs[[i]], extractMarker(batch.out[[i]]$exprs, m))
        }
    }

    # Checking for correct behaviour with NAs.
    batch.weights[[1]][1:5] <- NA
    batch.weights[[2]][6:10] <- NA
    batch.weights[[3]][c(1,10)] <- NA
    FUN <- cydar:::.pullOutMarkers(batch.out, batch.weights)

    ref.weights <- batch.weights
    for (i in seq_along(ref.weights)) {
        cur.weights <- batch.weights[[i]]
        keep <- !is.na(cur.weights)
        ncells <- vapply(batch.out[[i]]$exprs[keep], nrow, 0L)
        ref.weights[[i]] <- rep(cur.weights[keep], ncells)
    }

    for (m in marker.names) {
        test <- FUN(m)
        expect_identical(ref.weights, test$weights)
        for (i in seq_along(batch.out)) {
            cur.weights <- batch.weights[[i]]
            keep <- !is.na(cur.weights)
            expect_identical(test$exprs[[i]], extractMarker(batch.out[[i]]$exprs[keep], m))
        }
    }
})

###########################################################

test_that("Range-based normalization functions are correct", {
    weighted.quantile <- function(obs, wts, p) {
        o <- order(obs)
        wts <- wts[o]
        obs <- obs[o]
        mid.cum.weight <- cumsum(wts) - wts/2
        total.weight <- sum(wts)
        approx(mid.cum.weight/total.weight, obs, xout=p, rule=2)$y
    }

    # Putting together observations.
    nbatches <- length(batch.out)
    M <- 2
    all.obs <- vector("list", nbatches)
    for (b in seq_len(nbatches)) { 
        cur.out <- batch.out[[b]]
        nsamples <- length(cur.out$exprs)
        
        cur.obs <- vector("list", nsamples)
        for (s in seq_len(nsamples)) { 
            cur.obs[[s]] <- cur.out$exprs[[s]][,M]
        }
        all.obs[[b]] <- unlist(cur.obs)
    }

    Q <- mapply(weighted.quantile, all.obs, test.wts, MoreArgs=list(p=c(0.01, 0.99)))
    all.FUN <- cydar:::.rescaleDistr(all.obs, test.wts, target=NULL, p=0.01)
    for (b in seq_len(nbatches)) { 
        expect_equal(rowMeans(Q), all.FUN[[b]](Q[,b]))
    }

    Q <- mapply(weighted.quantile, all.obs, test.wts, MoreArgs=list(p=c(0.05, 0.95)))
    all.FUN <- cydar:::.rescaleDistr(all.obs, test.wts, target=NULL, p=0.05)
    for (b in seq_len(nbatches)) { 
        expect_equal(rowMeans(Q), all.FUN[[b]](Q[,b]))
    }

    # With target set.
    Q <- mapply(weighted.quantile, all.obs, test.wts, MoreArgs=list(p=c(0.05, 0.95)))
    all.FUN <- cydar:::.rescaleDistr(all.obs, test.wts, target=2, p=0.05)
    for (b in seq_len(nbatches)) { 
        expect_equal(Q[,2], all.FUN[[b]](Q[,b]))
    }
})

test_that("Overall normalization function does its job", {
    get.meanvar <- function(stuff) {
        collected <- list()
        for (b in seq_along(stuff)) {
            collected[[b]] <- sapply(stuff[[b]], colMeans)
        }
        collected <- do.call(rbind, collected)
        apply(collected, 2, var)
    }

    # Computing variance in the means for all batches.
    oldV <- get.meanvar(all.x)
    post.xr <- normalizeBatch(all.x, batch.comp)
    expect_identical(names(post.xr), names(all.x))
    expect_identical(lapply(post.xr, names), lapply(all.x, names))
    
    post.xrU <- unlist(post.xr, recursive=FALSE)
    expect_identical(lapply(unlist(all.x, recursive=FALSE), colnames), lapply(post.xrU, colnames))

    newV <- get.meanvar(post.xr)
    expect_true(all(oldV > newV))

    # Repeating with warping, using only the first two markers for speed.
    set.seed(9999)
    post.xw <- normalizeBatch(all.x, batch.comp, markers=c("X1", "X2"), mode="warp")
    newV <- get.meanvar(post.xw)
    expect_true(all(oldV[1:2] > newV))

    # Checking that setting markers does the same thing.
    chosen <- c("X2", "X7")
    post.xm <- normalizeBatch(all.x, batch.comp, markers=chosen)
    expect_identical(names(post.xm), names(all.x))
    expect_identical(lapply(post.xm, names), lapply(all.x, names))

    post.xmU <- unlist(post.xm, recursive=FALSE)
    expect_true(all(unlist(lapply(post.xmU, FUN=function(y) { identical(colnames(y), chosen) }))))
    for (i in seq_along(post.xmU)) { 
        expect_equal(post.xmU[[i]][,chosen], post.xrU[[i]][,chosen])
    }
})

