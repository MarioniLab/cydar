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
# Note that downsampling is done _within_ each batch, hence the use of 'cell.id' via .raw_cell_index().
# We then have to use 'ordering' as 'i' is relative to the reordered indices, and we need it on the merged indices.
{
    ordering <- .raw_precomputed(x)$order
    i <- ((.raw_cell_index(x)[ordering] - 1L) %% as.integer(downsample))==0L
    ordering[i]
}
