#' @export
#' @importFrom graphics plot
plotCellLogFC <- function(x, y, logFC, max.logFC=NULL, zero.col=0.8, length.out=100, pch=16, ...) 
# Visualizes the cells with appropriate coloration, using a PCA plot by default.
#
# written by Aaron Lun
# created 20 April 2016
{
    if (is.null(max.logFC)) {
        max.logFC <- max(abs(logFC))
    } else {
        logFC[logFC < -max.logFC] <- -max.logFC
        logFC[logFC > max.logFC] <- max.logFC
    }

    # Linear interpolator of colours.
    logFC.col <- .interpolate_cols(logFC, max.logFC, left=c(0,0,1), middle=rep(zero.col, 3), right=c(1,0,0))

    if (pch %in% 21:25) {           
        plot(x, y, pch=pch, bg=logFC.col, ...)
    } else {
        plot(x, y, pch=pch, col=logFC.col, ...)
    }
                          
    values <- seq(from=-max.logFC, to=max.logFC, length.out=length.out)
    stored <- .interpolate_cols(values, max.logFC, left=c(0,0,1), middle=rep(zero.col, 3), right=c(1,0,0))
    names(stored) <- values
    return(invisible(stored))
}

#' @importFrom grDevices rgb
.interpolate_cols <- function(val, max.val, left, middle, right) 
# Need to account for the fact that RGB colors are squared.
{
    prop <- val/max.val
    prop[prop < -1] <- -1
    prop[prop > 1] <- 1

    pos <- prop > 0 
    pos.output <- neg.output <- vector("list", 3)

    for (i in seq_along(pos.output)) {
        true.right <- right[i]^2
        true.left <- left[i]^2
        true.middle <- middle[i]^2

        pgrad <- true.right - true.middle
        pos.output[[i]] <- sqrt(pgrad * prop[pos] + true.middle)

        ngrad <- true.left - true.middle
        neg.output[[i]] <- sqrt(ngrad * -prop[!pos] + true.middle)
    }

    all.cols <- character(length(pos))
    all.cols[pos] <- rgb(pos.output[[1]], pos.output[[2]], pos.output[[3]])
    all.cols[!pos] <- rgb(neg.output[[1]], neg.output[[2]], neg.output[[3]])
    return(all.cols)
}

#' @export
#' @importFrom graphics plot
#' @importFrom viridis viridis
plotCellIntensity <- function(x, y, intensity, irange=NULL, length.out=100, pch=16, ...) 
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

    all.cols <- viridis(n=length.out)
    mdpts <- seq(from=irange[1], to=irange[2], length.out=length.out)

    # Note that 'mdpts' represents the midpoints of each category, 
    # so we have to add half the category width to get the actual boundary.
    actual.threshold <- mdpts + (mdpts[2] - mdpts[1])/2
    ix <- pmin(length.out, findInterval(intensity, actual.threshold) + 1)
    cur.cols <- all.cols[ix]

    if (pch %in% 21:25) {           
        plot(x, y, pch=pch, bg=cur.cols, ...)
    } else {
        plot(x, y, pch=pch, col=cur.cols, ...)
    }
                          
    names(all.cols) <- mdpts 
    return(invisible(all.cols))
}
