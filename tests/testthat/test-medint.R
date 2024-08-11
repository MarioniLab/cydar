# Testing the median intensity function.
# require(cydar); require(testthat); source("test-medint.R")

set.seed(100)
nmarkers <- 10

# Setup.
ncells1 <- 1001
all.values1 <- matrix(rnorm(ncells1*nmarkers, sd=1), nrow=ncells1, ncol=nmarkers)
colnames(all.values1) <- paste0("X", seq_len(nmarkers))

ncells2 <- 2001
all.values2 <- matrix(rnorm(ncells2*nmarkers, sd=1), nrow=ncells2, ncol=nmarkers)
colnames(all.values2) <- colnames(all.values1)
fs <- list(A=all.values1, B=all.values2)

test_that("medIntensities works as expected", {
    to.use <- rep(c(TRUE, FALSE), c(8, 2))
    suppressWarnings(cd <- prepareCellData(fs, markers=to.use))
    out <- countCells(cd, downsample=1L, filter=1L)
    rcnt <- medIntensities(out, markers=!to.use)

    # Median intensities are computed correctly.
    combined <- rbind(all.values1[,to.use], all.values2[,to.use])
    ref.groups <- BiocNeighbors::findNeighbors(combined, threshold=int_metadata(out)$cydar$tol * sqrt(sum(to.use)), get.distance=FALSE)$index
    sample.ids <- rep(c("A", "B"), c(ncells1, ncells2))

    for (u in which(!to.use)) { 
        cur.assay <- assay(rcnt, paste0("med.", markernames(rcnt, mode="all")[u]))
        ref.assay <- matrix(NA_real_, length(ref.groups), ncol(rcnt))
        colnames(ref.assay) <- colnames(rcnt)

        all.int <- c(all.values1[,u], all.values2[,u])
        for (r in seq_along(ref.groups)) { 
            cur.group <- union(r, ref.groups[[r]]) # adding back self.
            by.sample <- split(cur.group, sample.ids[cur.group])

            for (s in names(by.sample)) { 
                ref.assay[r,s] <- median(all.int[by.sample[[s]]])
            }
        }
        expect_equal(ref.assay, cur.assay)
    } 

    # Other fields remain unchanged.
    expect_identical(assay(out, "counts"), assay(rcnt, "counts"))
    expect_identical(cellAssignments(out), cellAssignments(rcnt))
    expect_identical(intensities(out), intensities(rcnt))

    # Responds to downsampling correctly.
    alt <- countCells(cd, downsample=6L, filter=1L)
    rcnt2 <- medIntensities(alt, markers=!to.use)
    expect_equal(rcnt2, rcnt[c(seq_len(ncells1), seq_len(ncells2)) %% 6L == 1L,])

    # Handles alternative inputs.
    rcnt3 <- medIntensities(out, markers=c("X9", "X10"))
    expect_equal(rcnt, rcnt3)

    # Handles empty inputs.
    rcnt4 <- medIntensities(out, markers=character(0))
    expect_equal(out, rcnt4)
})
