#' @export
#' @importFrom BiocNeighbors findKNN
neighborDistances <- function(x, neighbors=50, downsample=50, as.tol=TRUE)
# Calculates the 'tol' required to capture a certain number of neighbors.
#
# written by Aaron Lun
# created 7 July 2016
{
    pre <- x$precomputed
    to.check <- .downsample0(x$cell.id, downsample)

    # Computing distances.
    distances <- findKNN(BNINDEX=pre, k=neighbors, get.index=FALSE, subset=to.check, raw.index=TRUE)$distance

    # Converting to tolerance values, if so desired.
    if (as.tol) {
        distances <- distances/sqrt(ncol(pre))
    }
    distances
}
