#' @export
#' @importFrom BiocNeighbors findNeighbors buildIndex
#' @importFrom methods is
findFirstSphere <- function(x, pvalues, threshold=1, block=NULL, num.threads=1)
# Returns a logical vector indicating which hyperspheres are redundant
# within the specified distance threshold.
#
# written by Aaron Lun
# created 31 October 2016
{
    if (length(pvalues) != nrow(x)) {
        stop("length of 'pvalues' must equal number of cells in 'x'")
    }

    if (is(x, "CyData")) {
        .check_cell_data(x)
        x <- .raw_intensities(x) # hyperspheres are in the rows.
    }

    if (!is.null(block)) {
        # Identifying unique elements within each block.
        if (length(block) != nrow(x)) {
            stop("length of 'block' must equal number of rows in 'x'")
        }
        by.block <- split(seq_along(block), block)
        total.out <- logical(length(block))
        for (b in by.block) {
            total.out[b] <- Recall(x[b,,drop=FALSE], pvalues[b], threshold=threshold, block=NULL)
        }
        return(total.out)
    }

    # Using findNeighbors to screen out candidates based on the hypersphere
    # that encloses the hypercube that defines the minimum distance boundary.
    pre <- buildIndex(x)
    MULT <- max(1, sqrt(ncol(x)))
    potential <- findNeighbors(pre, threshold=threshold * MULT, get.distance=FALSE, num.theads=num.threads)$index

    # Transposing for more efficient.
    drop_redundant(t(x), order(pvalues) - 1L, potential, threshold)
}

