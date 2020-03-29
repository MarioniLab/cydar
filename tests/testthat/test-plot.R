# Testing the plotting functions.
# require(cydar); require(testthat); source("test-plot.R")

set.seed(200)
x <- runif(1000)
y <- runif(1000)

test_that("plotSphereLogFC works correctly", {
    lfc <- rnorm(1000)
    
    out <- plotSphereLogFC(x, y, lfc)
    maxchange <- max(abs(lfc))
    expect_equal(maxchange, -as.double(head(names(out), 1)))
    expect_equal(maxchange, as.double(tail(names(out), 1)))
    
    out <- plotSphereLogFC(x, y, lfc, max.logFC=2)
    maxchange <- 2
    expect_equal(maxchange, -as.double(head(names(out), 1)))
    expect_equal(maxchange, as.double(tail(names(out), 1)))
})

test_that("plotSphereIntensity works correctly", {
    intensities <- rgamma(1000, 2, 2)
    out <- plotSphereIntensity(x, y, intensities)
    expect_equal(min(intensities), as.double(head(names(out), 1)))
    expect_equal(max(intensities), as.double(tail(names(out), 1)))
    
    out <- plotSphereIntensity(x, y, intensities, irange=c(1, 4))
    expect_equal(1, as.double(head(names(out), 1)))
    expect_equal(4, as.double(tail(names(out), 1)))
})

test_that("intensityRanges works correctly", {
    stuff <- list(A=matrix(rgamma(10000, 2, 2), ncol=20))
    colnames(stuff$A) <- paste0("X", seq_len(ncol(stuff$A)))
    cd <- prepareCellData(stuff)
    cn <- countCells(cd)
    
    x <- intensityRanges(cn, p=0.01)
    ref <- apply(stuff$A, 2, quantile, p=c(0.01, 0.99))
    rownames(ref) <- c("min", "max")
    expect_equal(x, ref)
    
    x <- intensityRanges(cn, p=0.05)
    ref <- apply(stuff$A, 2, quantile, p=c(0.05, 0.95))
    rownames(ref) <- c("min", "max")
    expect_equal(x, ref)

    # Works when there are only a subset of markers in use.
    cd <- prepareCellData(stuff, markers=c("X1", "X10"))
    alt <- intensityRanges(cn, p=0.05)
    expect_identical(x[,colnames(alt)], alt)
    expect_identical(colnames(alt), markernames(cn, mode="all"))
})
