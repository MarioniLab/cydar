# Tests expandRadius().
# library(cydar); library(testthat); source("test-expand.R")

set.seed(100000)
nmarkers <- 10
tol <- 0.5

# Setup.
fs <- list()
for (x in 1:10) { 
    ncells <- round(runif(1, 500, 2000))
    all.values <- matrix(rnorm(ncells*nmarkers, sd=1), nrow=ncells, ncol=nmarkers)
    colnames(all.values) <- paste0("Y", seq_len(nmarkers))
    fs[[paste0("x", x)]] <- all.values
}

cd <- prepareCellData(fs)

test_that("expandRadius works correctly", {
    all.means <- do.call(rbind, lapply(fs, colMeans))

    # No design matrix.
    out <- apply(all.means, 2, FUN=function(y) {
        fit <- lm(y ~ 1)
        summary(fit)$sigma^2
    })
    ex <- expandRadius(cd, tol=0.5)
    expect_equal(ex, sqrt(2 * mean(out) + 0.5^2))
    ex <- expandRadius(cd, tol=0)
    expect_equal(ex, sqrt(2 * mean(out)))

    # With a design matrix.
    g <- gl(5, 2)
    out <- apply(all.means, 2, FUN=function(y) {
        fit <- lm(y ~ g)
        summary(fit)$sigma^2
    })
    ex <- expandRadius(cd, design=model.matrix(~g), tol=0)
    expect_equal(ex, sqrt(2 * mean(out)))

    # With another design matrix.
    b <- factor(rep(1:5, 2))
    out <- apply(all.means, 2, FUN=function(y) {
        fit <- lm(y ~ g + b)
        summary(fit)$sigma^2
    })
    ex <- expandRadius(cd, design=model.matrix(~g + b), tol=0)
    expect_equal(ex, sqrt(2 * mean(out)))
})  
