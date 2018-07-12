#' @export
#' @importFrom kmknn findNeighbors precluster
#' @importFrom methods is
findFirstSphere <- function(x, pvalues, threshold=1, block=NULL)
# Returns a logical vector indicating which hyperspheres are redundant
# within the specified distance threshold.
#
# written by Aaron Lun
# created 31 October 2016
{
    if (length(pvalues)!=nrow(x)) {
        stop("length of 'pvalues' must equal number of rows in 'x'")
    }

    if (is(x, "CyData")) {
        .check_cell_data(x)
        x <- .raw_intensities(x)
    }

    if (!is.null(block)) {
        # Identifying unique elements within each block.
        if (length(block)!=nrow(x)) {
            stop("length of 'block' must equal number of rows in 'x'")
        }
        by.block <- split(seq_along(block), block)
        total.out <- logical(length(block))
        for (b in by.block) {
            total.out[b] <- Recall(x[b,,drop=FALSE], pvalues[b], threshold=threshold, block=NULL)
        }
        return(total.out)
    }

    # Using findNeighbors to screen out candidates based on the enclosing hypersphere.
    pre <- precluster(x)
    MULT <- max(1, sqrt(nrow(x)))
    potential <- findNeighbors(precomputed=pre, threshold=threshold * MULT, get.distance=FALSE, raw.index=TRUE)$index

    pvalues <- pvalues[pre$order]
    out <- .Call(cxx_drop_redundant, pre$data, order(pvalues) - 1L, potential, threshold)
    out[pre$order] <- out
    return(out)
}

