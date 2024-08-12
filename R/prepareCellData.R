#' Prepare mass cytometry data
#'
#' Convert single-cell marker intensities from a mass cytometry experiment into a format for efficient counting.
#' 
#' @param x A named list of numeric matrices, 
#' where each matrix corresponds to a sample and contains expression intensities for each cell (row) and each marker (column).
#'
#' Alternatively, a ncdfFlowSet object containing the same information.
#' @param markers A character vector containing the names of the markers to use in downstream analyses.
#' @param ... Additional arguments to pass to \code{\link{buildIndex}}.
#'
#' @details
#' This function constructs a neighbor search index from the marker intensities of each cell in one or more samples.
#' The precomputed index is used to speed up downstream nearest-neighbour searching,
#' avoiding redundant work from repeated calls to \code{\link{countCells}} (e.g., with different values of \code{tol}).
#' 
#' If \code{markers} is specified, only the selected markers will be used in the precomputation.
#' This restricts the markers that are used in downstream functions - 
#' namely, \code{\link{countCells}} and \code{\link{neighborDistances}}.
#' By default, \code{markers=NULL} which means that all supplied markers will be used.
#' 
#' Markers that are \emph{not} in \code{markers} will be ignored in distance calculations.
#' However, their intensities are still stored in the output object, for use in functions like \code{\link{medIntensities}}.
#' 
#' @return 
#' A \linkS4class{List} containing precomputed values for use in \code{\link{countCells}}.
#' This includes:
#' \itemize{
#' \item \code{precomputed}, a prebuilt index for the neighbor search.
#' \item \code{sample.id}, an integer vector specifying the sample of origin for each cell in \code{precomputed}.
#' \item \code{cell.id}, an integer vector specifying the original index in 
#' the corresponding sample of \code{x} for each cell in \code{precomputed}.
#' \item \code{unused}, a matrix of intensity values for markers \emph{not} in \code{markers}.
#' \item \code{colData}, a \linkS4class{DataFrame} containing per-sample statistics.
#' }
#'
#' @author Aaron Lun
#' 
#' @examples
#' ### Mocking up some data: ###
#' nmarkers <- 20
#' marker.names <- paste0("X", seq_len(nmarkers))
#' nsamples <- 8
#' sample.names <- paste0("Y", seq_len(nsamples))
#' 
#' x <- list()
#' for (i in sample.names) {
#'     ex <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
#'     colnames(ex) <- marker.names
#'     x[[i]] <- ex
#' }
#' 
#' ### Running the function: ###
#' cd <- prepareCellData(x)
#' cd
#'
#' @seealso
#' \code{\link{countCells}}, where the output of this function is used to obtain hypersphere counts.
#'
#' @export
#' @importFrom BiocNeighbors buildIndex
#' @importFrom methods as
#' @importFrom S4Vectors DataFrame List
prepareCellData <- function(x, markers=NULL, ...) {
    cell.data <- .pull_out_data(x)
    sample.names <-  cell.data$samples
    marker.names <- cell.data$markers
    exprs.list <- cell.data$exprs

    exprs <- do.call(rbind, exprs.list)
    ncells.per.sample <- vapply(exprs.list, FUN=nrow, FUN.VALUE=0L)
    sample.id <- rep(seq_along(exprs.list), ncells.per.sample)
    cell.id <- unlist(lapply(ncells.per.sample, seq_len), use.names=FALSE)

    used <- .chosen_markers(markers, marker.names)
    used.exprs <- t(exprs[,used,drop=FALSE])
  
    List(
        precomputed=buildIndex(used.exprs, transposed=TRUE, ...),
        sample.id=sample.id,
        cell.id=cell.id,
        used=used.exprs,
        unused=t(exprs[,!used,drop=FALSE]),
        colData=DataFrame(row.names=sample.names, totals=ncells.per.sample)
    )
}

#' @importFrom methods is
#' @importFrom Biobase sampleNames exprs
#' @importFrom BiocGenerics colnames
.pull_out_data <- function(x)
# Pulling out data so we don't have to rely on ncdfFlowSet input.
{
    if (is.list(x)) { 
        sample.names <- names(x)
        if (is.null(sample.names)) { 
            sample.names <- seq_along(x)
        }

        if (!length(x)) {
            stop("number of samples must be positive")
        } 

        marker.names <- colnames(x[[1]])
        for (i in sample.names) {
            x[[i]] <- as.matrix(x[[i]])
            tmp <- colnames(x[[i]])
            if (is.null(tmp)) { stop("column names must be labelled with marker identities"); }
            stopifnot(identical(tmp, marker.names))
        }
        expr.val <- unname(x)

    } else if (is(x, "ncdfFlowSet")) {
        sample.names <- sampleNames(x)
        marker.names <- colnames(x)
        by.sample <- seq_along(sample.names)
        expr.val <- lapply(by.sample, FUN=function(i) exprs(x[[i]]))
    } else {
        stop("'cell.data' must be a list or ncdfFlowSet object") 
    }

    list(samples=sample.names, markers=marker.names, exprs=expr.val)
}
