#' @export
#' @importFrom BiocNeighbors findDistance findNeighbors buildIndex
#' @importFrom methods is
spatialFDR <- function(x, pvalues, neighbors=50, bandwidth=NULL, num.threads=1)
# This controls the spatial FDR across a set of plot coordinates.
# Each point is weighted by the reciprocal of its density, based on the specified 'radius'.
# A frequency-weighted version of the BH method is then applied to the p-values.
#
# written by Aaron Lun
# created 23 May 2016
{
    if (is(x, "CyData")) {
        coords <- .raw_intensities(x)
    } else {
        coords <- x
    }

    if (length(pvalues)!=nrow(coords)) { 
        stop("coords 'nrow' and p-value vector length are not the same") 
    }

    # Discarding NA pvalues.
    haspval <- !is.na(pvalues)
    if (!all(haspval)) {
        coords <- coords[haspval,,drop=FALSE]
        pvalues <- pvalues[haspval]
    }

    pre <- buildIndex(coords)

    # Defining the bandwidth.        
    if (is.null(bandwidth)) { 
        neighbors <- as.integer(neighbors)
        if (neighbors == 0L) { 
            bandwidth <- 0 
        } else if (neighbors < 0L) { 
            stop("'neighbors' must be a non-negative integer") 
        } else { 
            # Figuring out the bandwidth for KDE, as the median of distances to the n-th neighbour.
            distances <- findDistance(pre, k=neighbors, num.threads=num.threads)
            bandwidth <- median(distances)
        }
    } else {
        bandwidth <- as.double(bandwidth)
    }

    if (is.na(bandwidth) || bandwidth <= 0) {
        warning("setting a non-positive bandwidth to a small offset")
        bandwidth <- 1e-8 
    }

    # Computing densities with a tricube kernel.
    dist2neighbors <- findNeighbors(pre, threshold=bandwidth, num.threads=num.threads, get.index=FALSE)$distance
    densities <- compute_density(dist2neighbors, bandwidth)
    w <- 1/densities

    # Computing a density-weighted q-value.
    o <- order(pvalues)
    pvalues <- pvalues[o]
    w <- w[o]

    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
    adjp <- pmin(adjp, 1)

    if (!all(haspval)) {
        refp <- rep(NA_real_, length(haspval))
        refp[haspval] <- adjp
        adjp <- refp
    }
    return(adjp)
}
