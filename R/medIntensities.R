#' @export
#' @importFrom flowCore markernames
#' @importFrom SummarizedExperiment assayNames assay<- assayNames<-
#' @importFrom S4Vectors metadata
medIntensities <- function(x, markers)
# Computes the median intensity for each sample, in each hypersphere,
# for each marker. The idea is to fit this to a linear model.
#
# written by Aaron Lun
# created 2 December 2016
{
    .check_cell_data(x)
    all.markers <- markernames(x, mode="all")
    selected <- all.markers[.chosen_markers(markers, all.markers)]
    
    all.leftovers <- markernames(x, mode="unused")
    if (!all(selected %in% all.leftovers)) {
        stop("markers used for counting cannot be used for computing median intensities")
    }
    keep <- all.leftovers %in% selected

    unused <- metadata(x)$cydar$unused[keep,,drop=FALSE]
    out <- .Call(cxx_compute_median_int, unused, ncol(x), .raw_sample_id(x), .raw_cellAssignments(x))
    
    chosen.leftovers <- all.leftovers[keep]
    old.names <- assayNames(x)
    for (j in seq_along(chosen.leftovers)) { 
        assay(x, length(old.names)+j) <- out[[j]]
    }
    assayNames(x) <- c(old.names, sprintf("med.%s", chosen.leftovers))
    return(x)
}

