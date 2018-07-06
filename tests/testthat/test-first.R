# This checks that the findFirstSphere function is operating properly. 
# require(cydar); require(testthat); source("test-first.R")

set.seed(400)
test_that("findFirstSphere is working correctly", {
    nhypers <- 1000
    nmarkers <- 10

    coords <- matrix(rnorm(nhypers*nmarkers, sd=1), ncol=nmarkers)
    pval <- rbeta(nhypers, 1, 10)
    tcoords <- t(coords)
    o <- order(pval) 
    block <- sample(3, nhypers, replace=TRUE)
    
    for (threshold in c(0.5, 1, 2)) {
        xkeep <- findFirstSphere(coords, pvalues=pval, threshold=threshold)

        # Checking that the exclusion is correct.
        ref <- logical(nhypers)
        locals <- integer(nhypers)
        for (j in seq_along(o)) {
            i <- o[j]
            if (j > 1) {
                if (!any(colSums(abs(tcoords[,ref,drop=FALSE] - tcoords[,i]) <= threshold)==nmarkers)) {
                    ref[i] <- TRUE
                }
            } else {
                ref[i] <- TRUE
            }
        }
        
        expect_identical(ref, xkeep)

        # Trying with a blocking factor.
        bkeep <- findFirstSphere(coords, pvalues=pval, threshold=threshold, block=block)
        all.out <- logical(nhypers)
        for (i in unique(block)) {
            chosen <- block==i
            all.out[chosen] <- findFirstSphere(coords[chosen,,drop=FALSE], pvalues=pval[chosen], threshold=threshold)
        }
        expect_identical(all.out, bkeep)
    }
})

test_that("findFirstSphere behaves sensibly with different inputs", {
    ncells <- 2000
    nmarkers <- 5

    # Trying CyData inputs.
    coords <- matrix(rnorm(ncells*nmarkers, sd=1), ncol=nmarkers)
    colnames(coords) <- paste0("X", seq_len(nmarkers))
    cd <- prepareCellData(list(A=coords))
    cnt <- countCells(cd, downsample=5, filter=1L)
    pval <- rbeta(nrow(cnt), 1, 10)

    for (threshold in c(0.5, 1, 2)) {
        out <- findFirstSphere(cnt, pval, threshold=threshold)
        ref <- findFirstSphere(intensities(cnt), pval, threshold=threshold)
        expect_equal(out, ref)
    }

    # Trying silly inputs.
    expect_error(findFirstSphere(cnt[0,], pval, threshold=threshold), "length")
    expect_identical(findFirstSphere(cnt[0,], pval[0], threshold=threshold), logical(0))
})
