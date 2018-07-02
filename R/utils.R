#' @importFrom flowCore markernames
.check_cell_data <- function(x) 
# Checks incoming cell data, that it was properly processed by prepareCellData.
{
    raw_int <- .raw_intensities(x)
    if (is.null(raw_int) || !is.matrix(raw_int) || !is.numeric(raw_int)) {
        stop("'intensities' should be a numeric matrix, run 'countCells()'")
    }

    if (ncol(raw_int)!=length(markernames(x))) {
        stop("dimensionality of intensity matrix is not equal to number of markers")
    }

    raw_ca <- .raw_cellAssignments(x)
    if (is.null(raw_ca) || !is.list(raw_ca)) {
        stop("'cellAssignments' should be a list, please run 'countCells()'")
    }

    raw_cc <- .raw_center_cell(x)
    if (is.null(raw_cc)) {
        stop("missing 'center.cell' in 'rowData(x)', please run 'countCells()'")        
    }

    invisible(NULL)
}

.chosen_markers <- function(chosen.markers, all.markers) 
# Identifying the markers that were chosen for use.
{
    if (!is.null(chosen.markers)) {
        used <- logical(length(all.markers))
        if (is.character(chosen.markers)) {
            chosen.markers <- match(chosen.markers, all.markers)
            if (any(is.na(chosen.markers))) {
                stop("specified 'markers' not in available set of markers")
            }
        }
        used[chosen.markers] <- TRUE
    } else {
        used <- rep(TRUE, length(all.markers))
    }
    return(used)
} 

.downsample <- function(x, downsample) 
# This determines which points should be used upon downsampling.
# There are three levels of indices that we need to consider:
#   i) the ordering of cells within each batch as originally supplied, stored as 'cell.id'.
#  ii) the ordering of cells in the merged matrix but _before_ reordering, stored as 'precomputed$order'.
# iii) the actual ordering of cells in the columns of the expression matrix returned by 'precluster()'.
# 
# Downsampling is done based on every 'nth' cell, using the ordering in (i).
# Note that 'cell.id' is already permuted to reflect the reordered matrix in (iii).
# Thus, 'i' refers to the columns of the reordered matrix that we would like to use for neighbour searching.
# However, findKNN() and friends require the indices of the matrix before reordering; hence the use of (ii).
{
    i <- ((.raw_cell_id(x) - 1L) %% as.integer(downsample))==0L
    .raw_precomputed(x)$order[i]
}
