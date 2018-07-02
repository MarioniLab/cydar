#' @export
#' @importFrom flowCore markernames
#' @importFrom SummarizedExperiment assayNames assay<- assayNames<-
medIntensities <- function(x)
# Computes the median intensity for each sample, in each hypersphere,
# for each marker. The idea is to fit this to a linear model.
#
# written by Aaron Lun
# created 2 December 2016
{
    .check_cell_data(x, check.clusters=FALSE)

    sample.id <- cellData(x)$sample.id - 1L # Get to zero indexing.
    ci <- .raw_cellIntensities(x)
    out <- .Call(cxx_compute_median_int, ci, ncol(x), sample.id, .raw_cellAssignments(x))
    
    used.markers <- markernames(x)
    old.names <- assayNames(x)
    for (j in seq_along(used.markers)) { 
        assay(x, length(old.names)+j) <- out[[j]]
    }
    assayNames(x) <- c(old.names, paste0("med.", used.markers))
    return(x)
}

