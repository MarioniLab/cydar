# This tests the distance calculating machinery.
# library(cydar); library(testthat); source("test-dist.R")

set.seed(600)
test_that("neighborDistances works correctly", {
    ncells <- 1000
    nmarkers <- 10
    coords <- matrix(rnorm(ncells*nmarkers, sd=1), nrow=ncells, ncol=nmarkers)
    colnames(coords) <- paste0("X", seq_len(nmarkers))
    suppressWarnings(cd <- prepareCellData(list(A=coords)))

    # Computing reference distances.
    refdist <- as.matrix(dist(t(cd$used)))
    nn <- 50L
    refbands <- apply(refdist, 1, function(x) { sort(x)[2:(nn+1)] })
    dimnames(refbands) <- NULL
    refbands <- t(refbands)

    stuff <- neighborDistances(cd, downsample=1, neighbors=nn, as.tol=FALSE)
    expect_equal(stuff, refbands)
        
    stuff2 <- neighborDistances(cd, downsample=1, neighbors=nn)
    expect_equal(stuff, stuff2*sqrt(nmarkers))

    # Increasing the downsampling.
    ds <- 10
    stuff3 <- neighborDistances(cd, downsample=ds, neighbors=nn)
    chosen <- seq_len(nrow(stuff2)) %% ds == 1L
    expect_equal(stuff2[chosen,,drop=FALSE], stuff3)
})
