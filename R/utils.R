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
# Find points should be used upon downsampling, based on original cell IDs.
# This ensures that the same points are chosen regardless of the samples used.
{
    which(.raw_cell_id(x) %% as.integer(downsample) == 1L)
}
