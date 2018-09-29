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
    expect_identical(markernames(cd, mode="unused"), setdiff(colnames(all.values1), chosen))

    cn <- countCells(cd)
    expect_identical(markernames(cn), chosen)
    expect_identical(markernames(cn, mode="all"), colnames(all.values1))
    expect_identical(markernames(cn, mode="unused"), setdiff(colnames(all.values1), chosen))
})

test_that("CyData cell information getters work as expected", {
    # Testing cellIntensities().
    cd <- prepareCellData(list(all.values1, all.values2))
    cn <- countCells(cd, filter=1)

    out <- cellIntensities(cd)
    expect_identical(out, cellIntensities(cn))
    expect_identical(out, cellIntensities(cn, mode="all"))
    expect_identical(dim(cellIntensities(cn, mode="unused")), c(0L, ncol(out)))

    # Testing with unused markers.
    chosen <- c("X2", "X7", "X8")
    cd2 <- prepareCellData(list(all.values1, all.values2), markers=chosen)

    out2 <- cellIntensities(cd2)
    o <- BiocNeighbors::KmknnIndex_clustered_order(int_metadata(cd2)$cydar$precomputed)
    ref <- t(rbind(all.values1, all.values2))[,o]
    expect_identical(out2, ref[chosen,])

    unused <- setdiff(rownames(ref), chosen)
    expect_identical(cellIntensities(cd2, mode="all"), ref[c(chosen, unused),])
    expect_identical(cellIntensities(cd2, mode="unused"), ref[unused,])

    # Testing cellInformation().
    info <- cellInformation(cd)
    expect_identical(nrow(info), ncol(out))

    o <- BiocNeighbors::KmknnIndex_clustered_order(int_metadata(cd)$cydar$precomputed)
    expect_identical(info$sample, rep(1:2, c(ncells1, ncells2))[o])
    expect_identical(info$row, c(seq_len(ncells1), seq_len(ncells2))[o])

    # Testing getCenterCell().
    centers <- getCenterCell(cn)
    index <- info$row[centers]
    expect_true(all(index%%10==1L))
    expect_error(getCenterCell(cd), "not available")
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
