#' @export
findFirstSphere <- function(x, pvalues, threshold=1, block=NULL, naive=FALSE)
# Returns a logical vector indicating which hyperspheres are redundant
# within the specified distance threshold.
#
# written by Aaron Lun
# created 31 October 2016
{
    if (length(pvalues)!=nrow(x)) {
        stop("length of 'pvalues' must equal number of rows in 'coords'")
    }

    if (!is.null(block)) {
        # Identifying unique elements within each block.
        if (length(block)!=nrow(coords)) {
            stop("length of 'block' must equal number of rows in 'coords'")
        }
        by.block <- split(seq_along(block), block)
        total.out <- logical(length(block))
        for (b in by.block) {
            total.out[b] <- Recall(coords[b,,drop=FALSE], pvalues[b], threshold=threshold, block=NULL, naive=naive)
        }
        return(total.out)
    }

    # Checking for non-redundancy.
    .check_cell_data(x)
    .Call(cxx_drop_redundant, order(pvalues) - 1L, .raw_center_cell(x), cellAssignments(x))
}

