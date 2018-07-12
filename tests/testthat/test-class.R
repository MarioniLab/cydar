# This tests the construction and function of the CyData class.
# require(testthat); require(cydar); source("test-class.R")

nmarkers <- 10
ncells1 <- 1001
all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
colnames(all.values1) <- paste0("X", seq_len(nmarkers))
ncells2 <- 2001
all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
colnames(all.values2) <- colnames(all.values1)

test_that("CyData getters work as expected", {
    # Some of the getters fail correctly.
    cd <- prepareCellData(list(all.values1, all.values2))
    expect_identical(markernames(cd), colnames(all.values1))
    expect_error(intensities(cd), "not available")
    expect_error(cellAssignments(cd), "not available")

    # All getters now work.
    cn <- countCells(cd)
    expect_identical(markernames(cn), colnames(all.values1))
    expect_type(intensities(cn), "double")
    expect_identical(dim(intensities(cn)), c(nrow(cn), as.integer(nmarkers)))
    expect_type(cellAssignments(cn), "list")

    # markernames() responds to subsetting.
    chosen <- c("X1", "X2")
    cd <- prepareCellData(list(all.values1, all.values2), markers=chosen)
    expect_identical(markernames(cd), chosen)
    expect_identical(markernames(cd, mode="all"), colnames(all.values1))

    cn <- countCells(cd)
    expect_identical(markernames(cn), chosen)
    expect_identical(markernames(cn, mode="all"), colnames(all.values1))
})

test_that("CyData column methods trigger warnings", {
    cd <- prepareCellData(list(all.values1, all.values2))
    cn <- countCells(cd)

    expect_warning(out <- cbind(cn, cn), "columns are not independent")
    expect_identical(assay(out), cbind(assay(cn), assay(cn)))

    expect_warning(out <- cn[,1:2], "columns are not independent")
    expect_identical(assay(out), assay(cn)[,1:2])

    blah <- cn
    expect_warning(blah[,2] <- blah[,1], "columns are not independent")
    XXX <- assay(cn, withDimnames=FALSE)
    XXX[,2] <- XXX[,1] 
    expect_identical(assay(blah, withDimnames=FALSE), XXX)
})
