#' @export
#' @importFrom graphics plot
#' @importFrom viridis viridis
plotSphereIntensity <- function(x, y, intensity, irange=NULL, col.range=viridis(100), pch=16, ...) 
# Visualizes the cells with appropriate coloration, using a PCA plot by default.
#
# written by Aaron Lun
# created 20 April 2016
{
    if (is.null(irange)) { 
        irange <- range(intensity)
    } else {
        intensity[intensity > irange[2]] <- irange[2]
        intensity[intensity < irange[1]] <- irange[1]
    }

    length.out <- length(col.range)
    mdpts <- seq(from=irange[1], to=irange[2], length.out=length.out)

    # Note that 'mdpts' represents the midpoints of each category, 
    # so we have to add half the category width to get the actual boundary.
    actual.threshold <- mdpts + (mdpts[2] - mdpts[1])/2
    ix <- pmin(length.out, findInterval(intensity, actual.threshold) + 1)
    cur.cols <- col.range[ix]

    if (pch %in% 21:25) {           
        plot(x, y, pch=pch, bg=cur.cols, ...)
    } else {
        plot(x, y, pch=pch, col=cur.cols, ...)
    }
                          
    names(col.range) <- mdpts 
    invisible(col.range)
}
