#' @export
#' @importFrom kmknn findNeighbors
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
countCells <- function(x, tol=0.5, BPPARAM=SerialParam(), downsample=10, filter=10, naive=FALSE)
# Extracts counts for each cell in a CyTOF experiment, based on the number of surrounding cells 
# from each sample, given a prepared set of expression values for all cells in each sample.
#
# written by Aaron Lun
# created 21 April 2016
# last modified 24 May 2017
{
    all.markers <- markernames(x)
    distance <- tol * sqrt(length(all.markers)) 
    if (distance <= 0) {
        warning("setting a non-positive distance to a small offset")
        distance <- 1e-8        
    }
    
    # Only collating hyperspheres around every '10th' cell, for speed.
    pre <- .raw_precomputed(x)
    downsample <- as.integer(downsample)
    chosen <- which(((pre$order - 1L) %% downsample) == 0L)

    ci <- findNeighbors(precomputed=pre, threshold=distance, BPPARAM=BPPARAM, 
        raw.index=TRUE, subset=chosen, get.distance=FALSE)$index
    
    # Computing the associated statistics.
    sample.id <- .raw_sample_id(x)
    out.stats <- .Call(cxx_compute_hyperstats, pre$data, ncol(x), sample.id - 1L, ci)
    out.counts <- out.stats[[1]]
    out.coords <- out.stats[[2]]

    # Creating a new object.
    output <- CyData(assays=list(counts=out.counts), colData=colData(x), metadata=metadata(x))
    intensities(output) <- out.coords
    cellAssignments(output) <- ci
    output$totals <- tabulate(sample.id, ncol(x))
    metadata(output)$cydar$tol <- tol
    return(output)
}
