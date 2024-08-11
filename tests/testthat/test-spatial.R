# Testing functions related to computing the spatial FDR.
# require(cydar); require(testthat); source("test-spatial.R")

set.seed(300)
nhypers <- 1000
nmarkers <- 10
coords <- matrix(rnorm(nhypers*nmarkers, sd=1), nrow=nhypers, ncol=nmarkers)
colnames(coords) <- paste0("X", seq_len(nmarkers))

test_that("Density calculation works correctly", {   
    d2n <- BiocNeighbors::findKNN(coords, k=50, get.index=FALSE)$distance
    bandwidth <- median(d2n[,50])
    refbands <- BiocNeighbors::findNeighbors(coords, threshold=bandwidth, get.index=FALSE)$distance
    densities <- cydar:::compute_density(refbands, bandwidth)

    refdist <- as.matrix(dist(coords))
    diag(refdist) <- 0
    weightmat <- 1 - (refdist/bandwidth)^3
    weightmat[weightmat < 0] <- 0
    refdens <- rowSums(weightmat^3)
    names(refdens) <- NULL
    expect_equal(refdens, densities)
})

test_that("spatialFDR works correctly overall", {
    pval <- rbeta(nhypers, 1, 10)
    nn <- 20
    x <- spatialFDR(coords, pval, neighbors=nn)

    # Testing that it responds to the bandwidth.
    d2n <- BiocNeighbors::findKNN(coords, k=nn, get.index=FALSE)$distance
    bandwidth <- median(d2n[,nn])
    x2 <- spatialFDR(coords, pval, bandwidth=bandwidth)
    expect_equal(x, x2)

    # Comparing to a straight-up implementation in R.
    refdist <- as.matrix(dist(coords))
    refbands <- apply(refdist, 1, function(x) { sort(x)[nn] })
    weightmat <- 1 - (refdist/bandwidth)^3
    weightmat[weightmat < 0] <- 0
    refdens <- rowSums(weightmat^3)

    o <- order(pval)
    pval <- pval[o]
    w <- 1/refdens
    w <- w[o]
    qval <- rev(cummin(rev(pval * sum(w)/cumsum(w))))
    qval[o] <- pmin(1, qval)
    names(qval) <- NULL
    expect_equal(x, qval)

    # Checking what happens when we add NAs in the pvalues.
    pval.2 <- pval
    discarded <- c(10, 15, 20, 50)
    pval.2[discarded] <- NA_real_
    nax <- spatialFDR(coords, pval.2, neighbors=nn)
    nax.2 <- spatialFDR(coords[-discarded,], pval[-discarded], neighbors=nn)
    expect_equal(nax[-discarded], nax.2)
    expect_true(all(is.na(nax[discarded])))
})

test_that("spatialFDR repsonds to other inputs", {
    pval <- rbeta(nhypers, 1, 10)

    # Handles CyData inputs.
    cd <- prepareCellData(list(A=coords))
    cn <- countCells(cd, downsample=1, filter=0)
    expect_equal(spatialFDR(cn, pval), spatialFDR(intensities(cn), pval))

    # Handles empty inputs.
    expect_warning(spatialFDR(coords[0,], pval[0]), "non-positive")
    expect_identical(spatialFDR(coords[0,], pval[0], bandwidth=1), numeric(0))
    
    # Handles mismatch.
    expect_error(spatialFDR(coords, pval[0]), "not the same")

    # Handles zero neighbors.
    expect_warning(x <- spatialFDR(coords, pval, neighbors=0), "bandwidth")
    expect_warning(x2 <- spatialFDR(coords, pval, bandwidth=0), "bandwidth")
    expect_equal(x, x2)
    expect_equal(x, p.adjust(pval, method="BH"))
})
