#' Count cells in high-dimensional space
#' 
#' Count the number of cells from each sample lying inside hyperspheres in high-dimensional space.
#' 
#' @param prepared A \linkS4class{List} object produced by \code{\link{prepareCellData}}.
#' @param tol A numeric scalar to be used as the scaling factor for the hypersphere radius.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization 
#' is to be performed in \code{\link{findNeighbors}}.
#' @param downsample An integer scalar specifying the frequency with which cells are sampled to form hyperspheres.
#' @param filter An integer scalar specifying the minimum count sum required to report a hypersphere.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @details
#' Consider that each cell defines a point in M-dimensional space (where M is the number of markers), based on its marker intensities.
#' This function constructs hyperspheres and counts the number of cells from each sample lying within each hypersphere.
#' In this manner, the distribution of cells across the space can be quantified.
#' For each hypersphere, cell counts for all samples are reported along with the median intensity across the counted cells for each marker.
#' 
#' Each hypersphere is centered on a cell to ensure that only occupied spaces are counted.
#' However, for high-density spaces, this can result in many redundant hyperspheres.
#' To reduce computational work, only a subset of cells are used to define hyperspheres.
#' The downsampling frequency is specified by \code{downsample}, e.g., only every 10th cell is used to make a hypersphere by default.
#' 
#' Each hypersphere also has a radius of \code{tol*sqrt(M)} (this relationship avoids loss of counts as M increases).
#' \code{tol} can be interpreted as the acceptable amount of deviation in the intensity of a single marker for a given subpopulation.
#' The default value of 0.5 means that, for any one marker, cells with +0.5 or -0.5 intensity will be counted into the same subpopulation.
#' This value is sensible as intensities are usually on a log-10 scale, such that a total of 10-fold variability in marker intensities is tolerated.
#' 
#' The coordinates are reported as (weighted) medians across all cells in each hypersphere.
#' Compared to the center, the median better reflects the location of the hypersphere if the cells are not distributed around the centre.
#' Each cell is weighted inversely proportional to the total number of cells in the corresponding sample.
#' This ensures that large samples do not dominate the median calculation.
#' 
#' All hyperspheres with count sums below \code{filter} are removed by default.
#' Such hyperspheres do not have enough counts (and thus, information) for downstream analyses.
#' Removing them reduces the amount of memory required to form the output matrix.
#' 
#' @return 
#' A \linkS4class{CyData} object containing the following information:
#' \describe{
#' \item{\code{counts}}{An integer matrix of counts for each hypersphere (row) and sample (column) in the \code{assays} slot.}
#' \item{\code{intensities}:}{A numeric matrix of median intensities for each hypersphere (row) and marker (column), 
#' accessible with the \code{\link{intensities}} function.}
#' \item{\code{cellAssignments}:}{A list of integer vectors specifying the cells contained within each hypersphere, 
#' accessible with the \code{\link{cellAssignments}} function.}
#' \item{\code{totals}:}{An integer vector specifying the total number of cells in each sample, 
#' stored as a field in the \code{colData}.}
#' \item{\code{center.cell}:}{An integer vector specifying the cell that is used as the centre of each hypersphere, 
#' accessible with the \code{\link{getCenterCell}} function.}
#' }
#' Contents of \code{prepared} are also stored in the \code{\link{int_metadata}} of the output object. 
#' 
#' @author
#' Aaron Lun
#' 
#' @references
#' Lun ATL, Richard AC, Marioni JC (2017). 
#' Testing for differential abundance in mass cytometry data. 
#' \emph{Nat. Methods}, 14, 7:707-709.
#' 
#' Samusik N, Good Z, Spitzer MH et al. (2016).
#' Automated mapping of phenotype space with single-cell data.
#' \emph{Nat. Methods} 13:493-496
#' 
#' @seealso
#' \code{\link{prepareCellData}}, to generate the input object \code{prepared}.
#'
#' @examples
#' example(prepareCellData, echo=FALSE)
#' downsample <- 10L
#' tol <- 0.5
#' 
#' cnt <- countCells(cd, filter=1, downsample=downsample, tol=tol)
#' cnt
#' 
#' @export
#' @importFrom BiocNeighbors findNeighbors
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom SingleCellExperiment int_elementMetadata int_metadata int_colData SingleCellExperiment
countCells <- function(prepared, tol=0.5, num.threads=1, BPPARAM=SerialParam(), downsample=10, filter=10) {
    bdx <- prepared$precomputed

    # Scaling to account for the number of used markers.
    distance <- tol * sqrt(nrow(prepared$used))
    if (distance <= 0) {
        warning("setting a non-positive distance to a small offset")
        distance <- 1e-8
    }
    
    # Only collating hyperspheres around every n-th cell, for speed.
    chosen <- .downsample0(prepared$cell.id, downsample)
    ci <- findNeighbors(
        bdx,
        threshold=distance,
        num.threads=num.threads,
        BPPARAM=BPPARAM,
        raw.index=TRUE,
        subset=chosen,
        get.distance=FALSE)$index

    # Adding self back into its list, sorting to reduce cache misses.
    for (i in seq_along(ci)) {
        ci[[i]] <- sort(c(chosen[i], ci[[i]]))
    }

    # Filtering out low-abundance hyperspheres to avoid creating large matrices.
    keep <- lengths(ci) >= filter
    ci <- ci[keep]
    chosen <- chosen[keep]
        
    # Computing the counts.
    sample.id <- prepared$sample.id - 1L
    nsamples <- nrow(prepared$colData)
    out.counts <- count_cells(ci, sample.id, nsamples)
    out.counts <- t(out.counts)

    # Computing the median intensities (transposing for column-major access).
    sample.weights <- 1/tabulate(prepared$sample.id, nbins=nsamples)

    med.used <- weighted_median_int(t(prepared$used), ci, sample.id, sample.weights)
    med.used <- t(med.used)
    med.unused <- weighted_median_int(t(prepared$unused), ci, sample.id, sample.weights)
    med.unused <- t(med.unused)

    # Creating a new object (creating a SCE first to circumvent the CyData validity check).
    output <- SingleCellExperiment(
        assays=list(counts=out.counts), 
        colData=prepared$colData
    )

    int_elementMetadata(output)$cydar <- DataFrame(
        intensities=I(med.used),
        cellAssignments=I(ci), 
        center.cell=chosen,
        unused=I(med.unused)
    )

    prepared$colData <- NULL
    prepared$precomputed <- NULL
    int_metadata(output)$cydar <- prepared
    int_metadata(output)$cydar$tol <- tol

    # Reordering for various historical reasons.
    output <- as(output, "CyData")
    output[order(.raw_sample_id(output)[chosen], .raw_cell_id(output)[chosen]),]
}
