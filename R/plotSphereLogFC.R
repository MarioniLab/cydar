#' @export
#' @importFrom graphics plot
#' @importFrom grDevices col2rgb
plotSphereLogFC <- function(x, y, logFC, max.logFC=NULL, zero.col="grey80", left.col="blue", right.col="red", length.out=100, pch=16, ...) 
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
    left.col <- col2rgb(left.col)[,1]/255
    zero.col <- col2rgb(zero.col)[,1]/255
    right.col <- col2rgb(right.col)[,1]/255

    logFC.col <- .interpolate_cols(logFC, max.logFC, left=left.col, middle=zero.col, right=right.col)

    if (pch %in% 21:25) {           
        plot(x, y, pch=pch, bg=logFC.col, ...)
    } else {
        plot(x, y, pch=pch, col=logFC.col, ...)
    }
                          
    values <- seq(from=-max.logFC, to=max.logFC, length.out=length.out)
    stored <- .interpolate_cols(values, max.logFC, left=left.col, middle=zero.col, right=right.col)
    names(stored) <- values
    return(invisible(stored))
}

#' @importFrom grDevices rgb
.interpolate_cols <- function(val, max.val, left, middle, right) 
# Need to account for the fact that RGB values are square-rooted 
# representations of the actual colour intensity.
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

