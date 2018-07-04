## This sets up the CyData reference class.
## This differs from a standard SE class only in that it raises warnings upon column subsetting and binding.

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setClass("CyData", contains="SummarizedExperiment")

#' @importFrom S4Vectors setValidity2
setValidity2("CyData", function(object) {
    msg <- character(0)
        
    raw_sid <- .raw_sample_id(object)
    if (is.null(raw_sid) || !all(raw_sid > 0L & raw_sid <= ncol(object))){ 
        msg <- c(msg, "missing or invalid sample IDs")
    }

    raw_cin <- .raw_cellIntensities(object)
    raw_markers <- markernames(object)
    if (!identical(length(raw_markers), nrow(raw_cin))) {
        msg <- c(msg, "mismatch in number of markers and dimensionality of cell intensity matrix")
    }

    raw_cid <- .raw_sample_id(object)
    if (is.null(raw_cid) || !all(raw_cid > 0L & raw_cid <= ncol(raw_cin))){ 
        msg <- c(msg, "missing or invalid cell IDs")
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
        .col_warning()
    }
    callNextMethod()
})

#' @export
setMethod("cbind", "CyData", function(..., deparse.level=1) {
    .col_warning()
    callNextMethod()
})

################################################
# Defining a stub constructor.

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
CyData <- function(...) {
    new("CyData", SummarizedExperiment(...))
}

scat <- function(fmt, vals=character(), exdent=2, ...) {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

#' @export
#' @importFrom methods show
#' @importFrom flowCore markernames
setMethod("show", signature("CyData"), function(object) {
    callNextMethod()
    scat("markers(%d): %s\n", markernames(object))
    cat(sprintf("cells: %i\n", ncol(.raw_precomputed(object)$data)))
})

################################################
## Defining some getters and setters for internal use.

#' @importFrom S4Vectors metadata
.raw_metadata <- function(x) metadata(x)$cydar

.raw_precomputed <- function(x) .raw_metadata(x)$precomputed

.raw_cellIntensities <- function(x) .raw_precomputed(x)$data

.raw_sample_id <- function(x) .raw_metadata(x)$sample.id

.raw_cell_id <- function(x) .raw_metadata(x)$cell.id

#' @importFrom SummarizedExperiment rowData
.raw_rowdata <- function(x) rowData(x)$cydar

.raw_center_cell <- function(x) .raw_rowdata(x)$center.cell

.raw_cellAssignments <- function(x) .raw_rowdata(x)$cellAssignments

.raw_intensities <- function(x) .raw_rowdata(x)$intensities

################################################
## Defining some getters for external use.

#' @export
setGeneric("cellAssignments", function(x) standardGeneric("cellAssignments"))

#' @export
setMethod("cellAssignments", "CyData", function(x) {
    val <- .raw_cellAssignments(x)
    if (is.null(val)) {
        stop("cell assignments are not available in 'x', run 'countCells()'")
    }
    names(val) <- rownames(x)
    return(val)
})

#' @export
setGeneric("intensities", function(x) standardGeneric("intensities"))

#' @export
#' @importFrom flowCore markernames
setMethod("intensities", "CyData", function(x) {
    val <- .raw_intensities(x)
    if (is.null(val)) {
        stop("intensities are not available in 'x', run 'countCells()'")
    }
    colnames(val) <- markernames(x)
    rownames(val) <- rownames(x)
    val
})

#' @export
#' @importFrom flowCore markernames
setMethod("markernames", "CyData", function(object, all=FALSE) {
    mdf <- metadata(object)$cydar$markers
    all.markers <- rownames(mdf)
    if (!all) {
        all.markers <- all.markers[mdf$used]
    }
    return(all.markers)
})
