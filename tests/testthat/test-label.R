# Tests labelSpheres.
# library(testthat); library(cydar); source("test-label.R")

set.seed(400)
test_that("labelSpheres is working correctly", {
    nhypers <- 1000
    nmarkers <- 10

    coords <- matrix(rnorm(nhypers*nmarkers, sd=1), ncol=nmarkers)
    all.labels <- character(nhypers)
    all.labels[sample(nhypers, length(LETTERS))] <- LETTERS
    
    re.label <- labelSpheres(coords, all.labels)
    out <- BiocNeighbors::queryKNN(query=coords, X=coords[all.labels!="",], k=2)
    expect_identical(re.label, all.labels[all.labels!=""][out$index[,1]])

    # Handles duplicate labels.
    coords <- matrix(rnorm(nhypers*nmarkers, sd=1), ncol=nmarkers)
    all.labels <- character(nhypers)
    all.labels[sample(nhypers, 50, replace=TRUE)] <- sample(LETTERS, 50, replace=TRUE)

    re.label <- labelSpheres(coords, all.labels)
    out <- BiocNeighbors::queryKNN(query=coords, X=coords[all.labels!="",], k=2)
    expect_identical(re.label, all.labels[all.labels!=""][out$index[,1]])

    # Handles solo labels.
    coords <- matrix(rnorm(nhypers*nmarkers, sd=1), ncol=nmarkers)
    all.labels <- character(nhypers)
    all.labels[5:10] <- "A"

    re.label <- labelSpheres(coords, all.labels)
    expect_identical(re.label, rep("A", nhypers))

    # Handles no labels.
    coords <- matrix(rnorm(nhypers*nmarkers, sd=1), ncol=nmarkers)
    all.labels <- character(nhypers)
    re.label <- labelSpheres(coords, all.labels)
    expect_identical(re.label, all.labels)
})

test_that("labelSpheres behaves sensibly with different inputs", {
    ncells <- 2000
    nmarkers <- 5

    # Trying CyData inputs.
    coords <- matrix(rnorm(ncells*nmarkers, sd=1), ncol=nmarkers)
    colnames(coords) <- paste0("X", seq_len(nmarkers))
    cd <- prepareCellData(list(A=coords))
    cnt <- countCells(cd, downsample=5, filter=1L)

    all.labels <- character(nrow(cnt))
    all.labels[sample(nrow(cnt), length(LETTERS))] <- LETTERS
    out <- labelSpheres(cnt, all.labels)
    ref <- labelSpheres(intensities(cnt), all.labels)
    expect_equal(out, ref)

    # Trying silly inputs.
    expect_error(labelSpheres(cnt[0,], all.labels), "not TRUE")
    expect_identical(labelSpheres(cnt[0,], all.labels[0]), character(0))
})
