# Testing the prepareCellData machinery.
# require(testthat); require(cydar); source("test-prep.R")

nmarkers <- 10
ncells1 <- 1001
ncells2 <- 2001
all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)

colnames(all.values1) <- paste0("X", seq_len(nmarkers))
colnames(all.values2) <- colnames(all.values1)

set.seed(90001)
test_that("prepareCellData works as expected", {
    out <- prepareCellData(list(X=all.values1, Y=all.values2))

    # Checking that the cells are correctly ordered.
    sid <- rep(1:2, c(ncells1, ncells2))
    cid <- c(seq_len(ncells1), seq_len(ncells2))
    expect_identical(sid, out$sample.id)
    expect_identical(cid, unname(out$cell.id))

    # Checking that the intensities are correctly reordered.
    expect_equivalent(rbind(all.values1, all.values2), t(out$used))
    expect_identical(dim(out$unused), c(0L, as.integer(ncells1+ncells2)))
})

set.seed(90001)
test_that("prepareCellData works with subsetted markers", {
    # Initial check.
    spec <- c(2,3,4)
    set.seed(100)
    out.sub <- prepareCellData(list(X=all.values1, Y=all.values2), markers=spec)
    set.seed(100)
    out.ref <- prepareCellData(list(X=all.values1[,spec], Y=all.values2[,spec]))

    tmp.sub <- out.sub
    tmp.ref <- out.ref
    tmp.sub$unused <- NULL
    tmp.ref$unused <- NULL
    expect_equal(tmp.sub, tmp.ref)

    # Ensuring that the unused fields are valid.
    expect_identical(out.sub$unused, t(rbind(all.values1, all.values2)[,-spec,drop=FALSE]))

    # Checking that we can subset by name as well.
    spec <- c("X2","X3","X4")
    set.seed(100)
    out.sub2 <- prepareCellData(list(X=all.values1, Y=all.values2), markers=spec)
    expect_equal(out.sub, out.sub2)
})

test_that("prepareCellData behaves with ncdfFlowSet inputs", {
    library(ncdfFlow)
    data(GvHD)
    fs <- GvHD[1:2]
    ncfs <- ncdfFlowSet(fs)

    set.seed(100)
    out <- prepareCellData(ncfs)

    thing <- lapply(seq_along(sampleNames(ncfs)), function(i) flowCore::exprs(ncfs[[i]]))
    names(thing) <- sampleNames(ncfs)

    set.seed(100)
    ref <- prepareCellData(thing)

    expect_equal(ref, out)

    # Checking that the cells are correctly ordered.
    ncells1 <- nrow(fs[[1]])
    ncells2 <- nrow(fs[[2]])
    sid <- rep(1:2, c(ncells1, ncells2))
    cid <- c(seq_len(ncells1), seq_len(ncells2))
    expect_identical(sid, out$sample.id)
    expect_identical(cid, out$cell.id)
})

set.seed(90002)
test_that("prepareCellData behaves with silly inputs", {
    # No cells, no problems.
    expect_error(out <- prepareCellData(list(X=all.values1[0,], Y=all.values2[0,])), NA)

    # Handles missing column names.        
    out <- prepareCellData(list(all.values1, all.values2))
    expect_identical(rownames(out$colData), as.character(1:2))

    # Throws upon no samples.
    expect_error(out <- prepareCellData(list()), "must be positive")

    # Throws when one matrix has mismatching markers.
    tmp <- all.values2
    colnames(tmp) <- sample(colnames(tmp))
    expect_error(prepareCellData(list(all.values1, tmp)), "not TRUE")
})
