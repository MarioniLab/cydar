#' @export
#' @importFrom flowCore markernames
#' @importFrom stats quantile
intensityRanges <- function(x, p=0.01)
# Computes the log-fold change in intensity.
#
# written by Aaron Lun
# created 22 April 2016
{
    used.markers <- markernames(x)
    used.ci <- .raw_cellIntensities(x)
    used.ranges <- vector("list", length(used.markers))
    for (m in seq_along(used.markers)) {
        used.ranges[[m]] <- quantile(used.ci[m,], p=c(p, 1-p))
    }
    names(used.ranges) <- used.markers

    unused.markers <- markernames(x, mode="unused")
    unused.ci <- .raw_unusedIntensities(x)
    unused.ranges <- vector("list", length(unused.markers))
    for (m in seq_along(unused.markers)) {
        unused.ranges[[m]] <- quantile(unused.ci[m,], p=c(p, 1-p))
    }
    names(unused.ranges) <- unused.markers

    all.ranges <- c(used.ranges, unused.ranges)
    output <- do.call(cbind, all.ranges)
    rownames(output) <- c("min", "max")
    output
}
