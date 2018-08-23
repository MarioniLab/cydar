## This sets up the CyData reference class.
## This differs from a standard SE class only in that it raises warnings upon column subsetting and binding.

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

    raw_cid <- .raw_sample_id(object)
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
# Defining show methods.

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

#' @importFrom SingleCellExperiment int_metadata
.raw_metadata <- function(x) int_metadata(x)$cydar

.raw_precomputed <- function(x) .raw_metadata(x)$precomputed

.raw_cellIntensities <- function(x) .raw_precomputed(x)$data

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
setMethod("markernames", "CyData", function(object, mode=c("used", "all", "unused")) {
    mode <- match.arg(mode)
    mdata <- .raw_metadata(object)$markers
    switch(mode, used=mdata$used,
        unused=mdata$unused,
        all=c(mdata$used, mdata$unused))
})

#' @export
#' @importFrom flowCore markernames
cellIntensities <- function(x, mode=c("used", "all", "unused")) {
    mode <- match.arg(mode)
    if (mode=="used") {
        val <- .raw_cellIntensities(x)
    } else if (mode=="all") {
        val <- rbind(.raw_cellIntensities(x), .raw_unusedIntensities(x))
    } else {
        val <- .raw_unusedIntensities(x)
    }
 
    rownames(val) <- markernames(x, mode=mode)
    return(val) 
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
