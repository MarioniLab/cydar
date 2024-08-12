#' Compute distances to neighbors
#'
#' Calculate the distances in high-dimensional space to the neighboring cells.
#' 
#' @param prepared A \linkS4class{List} object produced by \code{\link{prepareCellData}}.
#' @param neighbors An integer scalar specifying the number of neighbours.
#' @param downsample An integer scalar specifying the frequency with which cells are examined.
#' @param as.tol A logical scalar specifying if the distances should be reported as tolerance values.
#' @param num.threads Integer scalar specifying the number of threads to use.
#' 
#' @details
#' This function examines each cell at the specified downsampling frequency, and computes the Euclidean distances to its nearest neighbors.
#' If \code{as.tol=TRUE}, these distances are reported on the same scale as \code{tol} in \code{\link{countCells}}.
#' This allows users to choose a value for \code{tol} based on the output of this function.
#' Otherwise, the distances are reported without modification.
#' 
#' To visualize the distances/tolerances, one option is to use boxplots, as shown below.
#' Each boxplot represents the distribution of tolerances required for hyperspheres to contain a certain number of cells.
#' For example, assume that at least 20 cells in each hypersphere are needed to have sufficient power for hypothesis testing.
#' Now, consider all hyperspheres that are large enough to include the 19th nearest neighbour.
#' The average distance required to do so would be the median of the boxplot generated from the 19th column of the output.
#' 
#' % It's fair to say that most hyperspheres should contain around 20 cells.
#' % If we use the median distance to the 19th nearest neighbor, then - well - exactly 50% of hyperspheres would be above 20, and 50% would be below 20.
#' 
#' Another option is to examine the distribution of counts at a given tolerance/distance.
#' This is done by counting the number of hyperspheres with a particular number of nearest neighbors closer than the specified tolerance.
#' In this manner, the expected count distribution from setting a particular tolerance can be determined.
#' Note that the histogram is capped at \code{neighbors} to save time.
#' 
#' Note that, for each examined cell, its neighbors are identified from the full set of cells.
#' Downsampling only changes the rate at which cells are examined, for the sake of computational efficiency.
#' Neighbors are not identified from the downsampled set as this will inflate the reported distances.
#' 
#' @return
#' A numeric matrix of distances where each row corresponds to an examined cell and each column \code{i} corresponds to the \code{i}th closest neighbor.
#' 
#' @examples
#' example(prepareCellData, echo=FALSE)
#' 
#' distances <- neighborDistances(cd, as.tol=FALSE)
#' boxplot(distances, xlab="Neighbor", ylab="Distance")
#' 
#' # Making a plot to choose 'tol' in countCells().
#' distances <- neighborDistances(cd, as.tol=TRUE)
#' boxplot(distances, xlab="Neighbor", ylab="Tolerance")
#' 
#' required.count <- 20 # 20 cells per hypersphere 
#' med <- median(distances[,required.count-1]) 
#' segments(-10, med, required.count-1, col="dodgerblue")
#' segments(required.count-1, med, y1=0, col="dodgerblue")
#' 
#' # Examining the distribution of counts at a given 'tol' of 0.7.
#' # (Adding 1 to account for the cell at the centre of the hypersphere.)
#' counts <- rowSums(distances <= 0.7) + 1
#' hist(counts, xlab="Count per hypersphere")
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{prepareCellData}}, to generate the \code{prepared} object.
#'
#' \code{\link{countCells}}, where the choice of \code{tol} can be guided by the distance distributions.
#'
#' @export
#' @importFrom BiocNeighbors findKNN
neighborDistances <- function(prepared, neighbors=50, downsample=50, as.tol=TRUE, num.threads=1)
# Calculates the 'tol' required to capture a certain number of neighbors.
#
# written by Aaron Lun
# created 7 July 2016
{
    pre <- prepared$precomputed
    to.check <- .downsample0(prepared$cell.id, downsample)
    distances <- findKNN(pre, k=neighbors, get.index=FALSE, subset=to.check, num.threads=num.threads)$distance

    # Converting to tolerance values by accounting for the number of used markers.
    if (as.tol) {
        distances <- distances / sqrt(nrow(prepared$used))
    }
    distances
}
