#' CyData class and methods
#'
#' The CyData class is derived from the \linkS4class{SingleCellExperiment} class.
#' It is intended to store the cell counts for each group of cells (rows) for each sample (columns).
#' Groups are intended to be hyperspheres (see \code{\link{countCells}}) but could also be arbitrary clusters of cells.
#' It also stores the median intensities for each group and the identity of cells in the groups, parallel to the rows.
#' 
#' CyData objects should not be created directly by users.
#' The class has some strict validity conditions that are not easily satisfied by manual construction.
#' Users should rely on functions like \code{\link{prepareCellData}} and \code{\link{countCells}} to create the objects.
#' An overview of the CyData class and the available methods.
#' 
#' @section Getter functions for group-level data:
#' In the following code chunks, \code{x} or \code{object} are CyData objects.
#' \code{mode} is a string specifying the types of markers that should be returned;
#' this defaults to only those markers that are used in \code{\link{prepareCellData}},
#' but can also return the unused markers or all of them.
#' 
#' \itemize{
#' \item \code{intensities(x, mode=c("used", "all", "unused")} 
#' returns a numeric matrix of intensities for each group of cells (rows) and markers (columns).
#' Rows of the output matrix correspond to rows of \code{x}.
#' Values are returned for the markers specified by \code{mode} (see above).
#' \item \code{cellAssignments(x)} returns a list of integer vectors, 
#' where each vector corresponds to a row of \code{x} and contains the indices of the cells in that group.
#' Indices refer to columns of \code{cellIntensities(x)}.
#' \item \code{markernames(object, mode=c("used", "all", "unused"))}
#' returns a character vector of the marker names, depending on \code{mode} (see above).
#' \item \code{getCenterCell(x)} returns the index of the cell used at the center of each hypersphere.
#' }
#'
#' @section Getter functions for cell-level data:
#' In the following code chunks, \code{x} is a CyData object and \code{mode} is as previously described.
#'
#' \itemize{
#' \item \code{cellIntensities(x, mode=c("used", "all", "unused"))} 
#' returns a numeric matrix of intensities for each marker (row) and cell (column).
#' \item \code{cellInformation(x)} returns a \linkS4class{DataFrame} with one row per cell.
#' The \code{sample} field specifies the sample of origin for each cell,
#' while the \code{cell} field specifies the original row index of that cell in its original sample.
#' }
#'
#' @section Subsetting:
#' The subsetting and combining behaviour of CyData objects is mostly the same as that of \linkS4class{SingleCellExperiment} objects.
#' The only difference is that subsetting or combining CyData objects by column is not advisable.
#' Indeed, attempting to do so will result in a warning from the associated methods.
#' This is because the columns are usually not independent in contexts involving clustering cells across multiple samples.
#' If a sample is to be removed, it is more appropriate to do so in the function that generates the CyData object (usually \code{\link{prepareCellData}}).
#' 
#' @examples
#' example(countCells, echo=FALSE)
#' 
#' markernames(cnt)
#' head(intensities(cnt))
#' head(cellAssignments(cnt))
#' 
#' @author
#' Aaron Lun
#'
#' @name CyData
#' @aliases CyData-class
#' markernames
#' markernames,CyData-method
#' cellAssignments
#' cellInformation
#' intensities
#' cellIntensities
#' getCenterCell
#' [,CyData,ANY,ANY,ANY-method
#' [<-,CyData,ANY,ANY,CyData-method
#' cbind,CyData-method
#' show,CyData-method
NULL

#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setClass("CyData", contains="SingleCellExperiment")

#' @importFrom S4Vectors setValidity2
setValidity2("CyData", function(object) {
    msg <- character(0)
        
    raw_sid <- .raw_sample_id(object)
    Ncells <- length(raw_sid)
    if (is.null(raw_sid) || !all(raw_sid > 0L & raw_sid <= ncol(object))){ 
        msg <- c(msg, "missing or invalid sample IDs")
    }

    used_cin <- .raw_cellIntensities(object)
    used_markers <- markernames(object)
    if (!identical(length(used_markers), nrow(used_cin))) {
        msg <- c(msg, "mismatch in number of markers and dimensionality of cell intensity matrix")
    }

    unused_cin <- .raw_unusedIntensities(object)
    unused_markers <- markernames(object, mode="unused")
    if (!identical(length(unused_markers), nrow(unused_cin))) {
        msg <- c(msg, "mismatch in number of markers and dimensionality of unused intensity matrix")
    }

    if (is.null(object$totals) || any(object$totals < 0L)) {
        msg <- c(msg, "total number of cells per sample should be a positive integer")
    }

    raw_cid <- .raw_cell_id(object)
    if (is.null(raw_cid) || !all(raw_cid > 0L & raw_cid <= object$totals[raw_sid])){ 
        msg <- c(msg, "missing or invalid cell IDs")
    }

    if (Ncells != length(raw_cid) || Ncells!=ncol(used_cin) || Ncells!=ncol(unused_cin)) {
        msg <- c(msg, "number of cells is not consistent across fields of 'object'")
    }

    if (length(msg)) return(msg) 
    return(TRUE)
})

################################################
# Adding warning to column-level operations.

.col_warning <- function() {
    warning("columns are not independent in a CyData object.
Rather, rerun prepareCellData() with a subset of samples."); 
}

#' @export
setMethod("[", c("CyData", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(j)) { 
        .col_warning()
    }
    callNextMethod()
})

#' @export
setMethod("[<-", c("CyData", "ANY", "ANY", "CyData"), function(x, i, j, ..., value) {
    if (!missing(j)) {
        stop("column replacement is not supported for 'CyData' objects")
    }
    callNextMethod()
})

#' @export
setMethod("cbind", "CyData", function(..., deparse.level=1) {
    .col_warning()
    callNextMethod()
})

################################################
# Defining show methods.

#' @export
#' @importFrom methods show
#' @importFrom flowCore markernames
#' @importFrom S4Vectors coolcat
setMethod("show", signature("CyData"), function(object) {
    callNextMethod()
    coolcat("markers(%d): %s\n", markernames(object))
    cat(sprintf("cells: %i\n", ncol(.raw_cellIntensities(object))))
})

################################################
## Defining some getters and setters for internal use.

#' @importFrom SingleCellExperiment int_metadata
.raw_metadata <- function(x) int_metadata(x)$cydar

.raw_cellIntensities <- function(x) .raw_metadata(x)$used

.raw_unusedIntensities <- function(x) .raw_metadata(x)$unused

.raw_sample_id <- function(x) .raw_metadata(x)$sample.id

.raw_cell_id <- function(x) .raw_metadata(x)$cell.id

#' @importFrom SingleCellExperiment int_elementMetadata
.raw_rowdata <- function(x) int_elementMetadata(x)$cydar

.raw_cellAssignments <- function(x) .raw_rowdata(x)$cellAssignments

.raw_intensities <- function(x) .raw_rowdata(x)$intensities

################################################
## Defining some getters for external use.

#' @export
cellAssignments <- function(x) {
    val <- .raw_cellAssignments(x)
    if (is.null(val)) {
        stop("cell assignments are not available in 'x', run 'countCells()'")
    }
    names(val) <- rownames(x)
    val
}

.named_intensities <- function(x) {
    out <- .raw_intensities(x)
    colnames(out) <- .get_used_markers(x)
    out
}

.named_unused_intensities <- function(x) {
    out <- int_elementMetadata(x)$cydar$unused
    colnames(out) <- .get_unused_markers(x)
    out    
}

#' @export
intensities <- function(x, mode=c("used", "all", "unused")) {
    mode <- match.arg(mode)
    switch(mode, 
        used=.named_intensities(x),
        unused=.named_unused_intensities(x),
        all=cbind(.named_intensities(x), .named_unused_intensities(x))
    )
}

.get_used_markers <- function(x) {
    rownames(.raw_cellIntensities(x))
}

.get_unused_markers <- function(x) {
    rownames(.raw_unusedIntensities(x))
}

#' @export
#' @importFrom flowCore markernames
setMethod("markernames", "CyData", function(object, mode=c("used", "all", "unused")) {
    mode <- match.arg(mode)
    switch(mode, 
        used=.get_used_markers(object),
        unused=.get_unused_markers(object),
        all=c(.get_used_markers(object), .get_unused_markers(object))
    )
})

#' @export
cellIntensities <- function(x, mode=c("used", "all", "unused")) {
    mode <- match.arg(mode)
    switch(mode,
        used=.raw_cellIntensities(x),
        all=rbind(.raw_cellIntensities(x), .raw_unusedIntensities(x)),
        unused=.raw_unusedIntensities(x)
    )
}

#' @export
#' @importFrom S4Vectors DataFrame
cellInformation <- function(x) {
    DataFrame(sample=.raw_sample_id(x), row=.raw_cell_id(x))
}

#' @export
getCenterCell <- function(x) {
    out <- .raw_rowdata(x)$center.cell
    if (is.null(out)) {
        stop("center cell not available, run 'countCells()'")
    }
    names(out) <- rownames(x)
    out
}
