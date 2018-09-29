#' @export
#' @importFrom BiocNeighbors findNeighbors KmknnIndex_clustered_data
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom SingleCellExperiment int_elementMetadata int_metadata int_colData SingleCellExperiment
countCells <- function(x, tol=0.5, BPPARAM=SerialParam(), downsample=10, filter=10)
# Extracts counts for each cell in a CyTOF experiment, based on the number of surrounding cells 
# from each sample, given a prepared set of expression values for all cells in each sample.
#
# written by Aaron Lun
# created 21 April 2016
{
    all.markers <- markernames(x)
    distance <- tol * sqrt(length(all.markers)) 
    if (distance <= 0) {
        warning("setting a non-positive distance to a small offset")
        distance <- 1e-8
    }
    
    # Only collating hyperspheres around every '10th' cell, for speed.
    pre <- .raw_precomputed(x)
    chosen <- .downsample(x, downsample)
    ci <- findNeighbors(precomputed=pre, threshold=distance, BPPARAM=BPPARAM, 
        raw.index=TRUE, subset=chosen, get.distance=FALSE)$index

    # Filtering out low-abundance hyperspheres to avoid creating large matrices.
    keep <- lengths(ci) >= filter
    ci <- ci[keep]
    chosen <- chosen[keep]
        
    # Computing the associated statistics.
    sample.id <- .raw_sample_id(x)
    out.stats <- .Call(cxx_compute_hyperstats, KmknnIndex_clustered_data(pre), ncol(x), sample.id - 1L, ci)
    out.counts <- out.stats[[1]]
    out.coords <- out.stats[[2]]

    # Creating a new object (again, creating a SCE first to circumvent the CyData validity check).
    output <- SingleCellExperiment(assays=list(counts=out.counts), colData=colData(x), metadata=metadata(x))
    int_metadata(output) <- int_metadata(x)
    int_colData(output) <- int_colData(x)
    int_elementMetadata(output)$cydar <- DataFrame(intensities=I(out.coords), cellAssignments=I(ci), center.cell=chosen)
    int_metadata(output)$cydar$tol <- tol
    output <- as(output, "CyData")

    # Reordering for various historical reasons.
    output <- output[order(sample.id[chosen], .raw_cell_id(output)[chosen]),]
    return(output)
}
