# This tests the countCells functions.
# library(cydar); library(testthat); source("test-counts.R")

nmarkers <- 10
tol <- 0.5

ncells1 <- 1001
all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
colnames(all.values1) <- paste0("X", seq_len(nmarkers))

ncells2 <- 2001
all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
colnames(all.values2) <- colnames(all.values1)
fs <- list(A=all.values1, B=all.values2)

set.seed(10000)
test_that("countCells computes all hypersphere values correctly", {
    # Counting with specified parameters.
    suppressWarnings(cd <- prepareCellData(fs))
    cn <- countCells(cd, filter=0L, downsample=1L, tol=tol)

    # Checking that the center choices are correct.
    centers <- int_elementMetadata(cn)$cydar$center.cell
    sid <- cydar:::.raw_sample_id(cn)[centers]
    expect_identical(sid, rep(1:2, c(ncells1, ncells2)))
    cid <- cydar:::.raw_cell_id(cn)[centers]
    expect_identical(unname(cid), c(seq_len(ncells1), seq_len(ncells2)))

    # Checking that the cell IDs are correct.
    total <- rbind(all.values1, all.values2)
    threshold <- tol * sqrt(nmarkers)

    ref <- BiocNeighbors::queryNeighbors(X=total, query=total, threshold=threshold, get.distance=FALSE)$index
    ref <- lapply(ref, sort)
    obs <- cellAssignments(cn)
    expect_identical(ref, obs)

    # Checking that the counts are correct.
    collected <- matrix(0L, nrow(total), ncol(cn))
    for (r in seq_len(nrow(total))) {
        cell.indices <- ref[[r]]
        collected[r,] <- tabulate(as.integer(cell.indices > ncells1)+1L, nbin=2)
    }
    expect_identical(collected, assay(cn, withDimnames=FALSE))

    # Checking the median per hypersphere.
    collected.meds.upper <- collected.meds.lower <- list()
    for (i in seq_along(ref)) { 
        cell.indices <- ref[[i]]
        combined <- total[cell.indices,,drop=FALSE]
        w <- 1/c(ncells1, ncells2)[as.integer(cell.indices > ncells1) + 1L]

        # Need to calculate lower/upper bounds for the median, as numerical imprecision has big effects.
        cur.meds <- apply(combined, 2, FUN=function(x) { 
            o <- order(x)
            x <- x[o]
            w <- w[o]
            p <- cumsum(w)/sum(w)
            x[c(sum(p < 0.499999), sum(p < 0.500001))+1]
        })
        collected.meds.lower[[i]] <- cur.meds[1,]
        collected.meds.upper[[i]] <- cur.meds[2,]
    }

    med.coords <- intensities(cn)
    collected.meds.lower <- do.call(rbind, collected.meds.lower)
    collected.meds.upper <- do.call(rbind, collected.meds.upper)
    expect_identical(dim(collected.meds.lower), dim(med.coords))
    expect_identical(dim(collected.meds.upper), dim(med.coords))

    out.of.range <- collected.meds.lower > med.coords | collected.meds.upper < med.coords
    expect_true(!any(out.of.range))
})

set.seed(10001)
test_that("countCells yields reproducible hypersphere-level results", {
    suppressWarnings(cd <- prepareCellData(fs))
    ref <- countCells(cd, filter=0L, downsample=1L, tol=tol)

    suppressWarnings(cd <- prepareCellData(fs))
    ref2 <- countCells(cd, filter=0L, downsample=1L, tol=tol)

    expect_identical(assay(ref), assay(ref2))
    expect_identical(intensities(ref), intensities(ref2))
    expect_identical(colData(ref), colData(ref2))
})

set.seed(10001)
test_that("countCells handles marker selection", {
    fs0 <- lapply(fs, function(x) {
        colnames(x) <- paste0(colnames(x), "++")
        x
    })
    fs2 <- mapply(FUN=cbind, fs, fs0, SIMPLIFY=FALSE)

    suppressWarnings(cd <- prepareCellData(fs2, markers=colnames(fs[[1]])))
    cn <- countCells(cd, filter=0L, downsample=1L, tol=tol)

    suppressWarnings(cd <- prepareCellData(fs))
    ref <- countCells(cd, filter=0L, downsample=1L, tol=tol)
    suppressWarnings(cd <- prepareCellData(fs))
    ref2 <- countCells(cd, filter=0L, downsample=1L, tol=tol)
    
    expect_identical(assay(cn), assay(ref))
    expect_identical(intensities(cn), intensities(ref))
    expect_identical(unname(intensities(cn, "unused")), unname(intensities(cn)))
})

set.seed(10001)
test_that("countCells responds to filtering", {
    suppressWarnings(cd <- prepareCellData(fs))
    cn <- countCells(cd, filter=0L, downsample=1L, tol=tol)

    # Responds to filtering.
    cn2 <- countCells(cd, filter=2L, downsample=1L, tol=tol)
    expect_equal(cn2, cn[rowSums(assay(cn))>=2L,])

    cn5 <- countCells(cd, filter=5L, downsample=1L, tol=tol)
    expect_equal(cn5, cn[rowSums(assay(cn))>=5L,])
})

set.seed(10001)
test_that("countCells responds to downsampling", {
    suppressWarnings(cd <- prepareCellData(fs))
    cn <- countCells(cd, filter=0L, downsample=1L, tol=tol)

    # Handles further downsampling.
    cid <- cydar:::.raw_cell_id(cn)[int_elementMetadata(cn)$cydar$center.cell]

    cn3 <- countCells(cd, filter=0L, downsample=3L, tol=tol)
    expect_equal(cn3, cn[cid %% 3L == 1L,])

    cn10 <- countCells(cd, filter=0L, downsample=10L, tol=tol)
    expect_equal(cn10, cn[cid %% 10L == 1L,])

    cn16 <- countCells(cd, filter=0L, downsample=16L, tol=tol)
    expect_equal(cn16, cn[cid %% 16L == 1L,])
})

set.seed(10001)
test_that("countCells handles parallelization", {
    library(BiocParallel)
    suppressWarnings(cd <- prepareCellData(fs))
    cn <- countCells(cd, filter=0L, downsample=1L, tol=tol)
    
    cn.p <- countCells(cd, filter=0L, downsample=1L, tol=tol, BPPARAM=MulticoreParam(2))
    expect_equal(cn, cn.p)

    cn.p <- countCells(cd, filter=0L, downsample=1L, tol=tol, num.threads=3)
    expect_equal(cn, cn.p)
})

test_that("countCells behaves correctly with silly inputs", { 
    nmarkers <- 10
    tol <- 0.5

    # Setup.
    ncells1 <- 57
    all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
    colnames(all.values1) <- paste0("X", seq_len(nmarkers))
    
    ncells2 <- 731
    all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
    colnames(all.values2) <- colnames(all.values1)
    fs <- list(A=all.values1, B=all.values2)

    suppressWarnings(cd <- prepareCellData(fs))

    # Testing nonsensical inputs.
    suppressWarnings(out <- countCells(cd, filter=1L, tol=0))
    expect_true(all(rowSums(assay(out))==1L))
    out <- countCells(cd, filter=Inf)
    expect_identical(nrow(out), 0L)

    suppressWarnings(empty <- prepareCellData(lapply(fs, "[", i=0,, drop=FALSE)))
    out <- countCells(empty)
    expect_identical(nrow(out), 0L)
})
