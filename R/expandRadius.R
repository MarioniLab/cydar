#' Expand the hypersphere radius
#' 
#' Expands the hypersphere radius to account for intensity shifting between non-barcoded samples.
#' 
#' @param prepared A \linkS4class{List} object produced by \code{\link{prepareCellData}}.
#' @param design A numeric matrix specifying the experimental design.
#' @param tol A numeric scalar proportional to the hypersphere radius, see \code{\link{countCells}}.
#' 
#' @details
#' This function increases the hypersphere radius to account for random shifts in marker intensity between non-barcoded samples.
#' The required increase is estimated by taking the mean of all intensities for each marker in each sample;
#' computing the variance of the mean intensities across samples for each marker;
#' and taking the mean variance across all markers.
#' This is equivalent to the square of the extra distance between cells caused by intensity shifts between samples.
#' 
#' The estimated increase is added onto \code{tol}, and the returned value can be directly used in the \code{tol} argument of \code{\link{countCells}}.
#' This expands the hyperspheres to ensure that corresponding subpopulations in different samples are still counted together.
#' Otherwise, an intensity shift in one sample may move the cells in a subpopulation out of a hypersphere.
#' This will inflate the variability if it occurs between replicate samples, and introduce spurious differences if it occurs between samples in different conditions.
#' 
#' @return
#' A numeric scalar specifying a modified \code{tol} to use in \code{\link{countCells}}.
#' 
#' @examples
#' ### Mocking up some data: ###
#' nmarkers <- 20
#' marker.names <- paste0("X", seq_len(nmarkers))
#' nsamples <- 8
#' sample.names <- paste0("Y", seq_len(nsamples))
#' 
#' x <- list()
#' for (i in sample.names) {
#'     ex <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
#'     ex <- t(t(ex) + rnorm(nmarkers, 0, 0.25)) # Adding a shift per marker
#'     colnames(ex) <- marker.names
#'     x[[i]] <- ex
#' }
#' 
#' ### Running the function: ###
#' cd <- prepareCellData(x)
#' expandRadius(cd)
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{prepareCellData}}, to generate the required input.
#'
#' \code{\link{countCells}}, where the \code{tol} can be set to the output of this function.
#' 
#' @export
#' @importFrom stats lm.fit
expandRadius <- function(prepared, design=NULL, tol=0.5) {
    sample.id <- prepared$sample.id
    nsamples <- nrow(prepared$colData)

    # Computing mean intensities for all (used) markers in all samples.
    all.means <- vector("list", nsamples)
    for (s in seq_len(nsamples)) { 
        all.means[[s]] <- rowMeans(prepared$used[,sample.id==s,drop=FALSE])
    }
    all.means <- do.call(rbind, all.means)
    
    # Fitting a linear model to estimate variance of shift process per marker.
    if (is.null(design)) { 
        design <- matrix(1, nrow=nrow(all.means), ncol=1)
    }
    fit <- lm.fit(y=all.means, x=design)

    # Taking the mean variance across all markers, representing squared shift.
    shift2 <- mean(fit$effects[-fit$qr$pivot[seq_len(fit$qr$rank)],]^2) 
    
    # Adding to default tolerance to get per-marker radius.
    # This is done based on squared distance, hence the squaring and rooting.
    # Doubling the shift to account for distances between cells, not just distance from cell to mean.
    sqrt(2*shift2 + tol^2)
}
