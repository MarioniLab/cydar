#' @export
#' @importFrom BiocNeighbors findKmknn
neighborDistances <- function(x, neighbors=50, downsample=50, as.tol=TRUE)
# Calculates the 'tol' required to capture a certain number of neighbors.
#
# written by Aaron Lun
# created 7 July 2016
{
    pre <- .raw_precomputed(x)
    to.check <- .downsample(x, downsample)

    # Computing distances.
    distances <- findKmknn(precomputed=pre, k=neighbors, get.index=FALSE, subset=to.check, raw.index=TRUE)$distance

    # Converting to tolerance values, if so desired.
    if (as.tol) {
        distances <- distances/sqrt(length(markernames(x)))
    }
    return(distances)
}
