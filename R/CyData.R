## This sets up the CyData reference class.
## This differs from a standard SE class only in that it raises warnings upon column subsetting and binding.

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setClass("CyData", contains="SummarizedExperiment")

.col_warning <- function() {
    warning("columns are not independent in a CyData object.
Rather, rerun prepareCellData() with a subset of samples."); 
}

setMethod("[", c("CyData", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(j)) { 
        .col_warning()
    }
    callNextMethod()
})

setMethod("[<-", c("CyData", "ANY", "ANY", "CyData"), function(x, i, j, ..., value) {
    if (!missing(j)) { 
        .colWarning()
    }
    callNextMethod()
})

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
.raw_precomputed <- function(x) metadata(x)$cydar$precomputed

.raw_cellIntensities <- function(x) .raw_precomputed(x)$data

#' @importFrom S4Vectors metadata
.raw_sample_id <- function(x) metadata(x)$cydar$sample.id

#' @importFrom SummarizedExperiment rowData
.raw_cellAssignments <- function(x) rowData(x)$cellAssignments

#' @importFrom SummarizedExperiment rowData
.raw_intensities <- function(x) rowData(x)$intensities

setGeneric("cellAssignments<-", function(x, value) standardGeneric("cellAssignments<-"))

setReplaceMethod("cellAssignments", "CyData", function(x, value) {
    rowData(x)$cellAssignments <- value
    return(x)
})

setGeneric("intensities<-", function(x, value) standardGeneric("intensities<-"))

setReplaceMethod("intensities", "CyData", function(x, value) {
    rowData(x)$intensities <- value
    return(x)
})

################################################
## Defining some getters for external use.

#' @export
setGeneric("cellAssignments", function(x) standardGeneric("cellAssignments"))

#' @export
setMethod("cellAssignments", "CyData", function(x) {
    out <- .raw_cellAssignments(x)
    names(out) <- rownames(x)
    return(out)
})

#' @export
setGeneric("intensities", function(x) standardGeneric("intensities"))

#' @export
#' @importFrom flowCore markernames
setMethod("intensities", "CyData", function(x) {
    val <- .raw_intensities(x)
    colnames(val) <- markernames(x)
    rownames(val) <- rownames(x)
    val
})

#' @export
#' @importFrom flowCore markernames
setMethod("markernames", "CyData", function(object) {
    return(metadata(object)$cydar$markers)
})


