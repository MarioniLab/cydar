#' @export
#' @importFrom kmknn precluster
prepareCellData <- function(x, markers=NULL, ...) 
# Converts it into a format more suitable for high-speed analysis.
# Also does k-means clustering to generate the necessary clusters.
#
# written by Aaron Lun
# created 14 August 2016
{
    cell.data <- .pull_out_data(x)
    sample.names <-  cell.data$samples
    marker.names <- cell.data$markers
    exprs.list <- cell.data$exprs

    exprs <- do.call(rbind, exprs.list)
    sample.id <- rep(seq_along(exprs.list), sapply(exprs.list, nrow))
    cell.id <- unlist(lapply(exprs.list, FUN=function(exprs) seq_len(ncol(exprs)) ))

    # Picking markers to use.
    used <- .chosen_markers(markers, marker.names)
    marker.names <- marker.names[used]
    reorg <- precluster(exprs[,used], ...)
  
    # Collating the output. 
    output <- CyData(colData=DataFrame(row.names=sample.names))
    metadata(output)$cydar <- list(
        precomputed=reorg,
        markers=marker.names[used],
        sample.id=sample.id[reorg$order],
        sample.id=cell.id[reorg$order]
    )
    output
}

#' @importFrom methods is
#' @importClassesFrom ncdfFlow ncdfFlowSet
#' @importFrom Biobase sampleNames 
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore exprs
.pull_out_data <- function(x)
# Pulling out data so we don't have to rely on ncdfFlowSet input.
{
    if (is.list(x)) { 
        sample.names <- names(x)
        if (is.null(sample.names)) { stop("list must be named by sample") }
        marker.names <- colnames(x[[1]])
        for (i in sample.names) {
            x[[i]] <- as.matrix(x[[i]])
            tmp <- colnames(x[[i]])
            if (is.null(tmp)) { stop("column names must be labelled with marker identities"); }
            stopifnot(identical(tmp, marker.names))
        }
        expr.val <- x
    } else if (is(x, "ncdfFlowSet")) { 
        sample.names <- Biobase::sampleNames(x)
        marker.names <- BiocGenerics::colnames(x)
        by.sample <- seq_along(sample.names)
        expr.val <- lapply(by.sample, FUN=function(i) flowCore::exprs(x[[i]]))
    } else {
        stop("'cell.data' must be a list or ncdfFlowSet object") 
    }
    return(list(samples=sample.names, markers=marker.names, exprs=expr.val))
}
