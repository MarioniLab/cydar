#' @export
#' @importFrom stats lm.fit
expandRadius <- function(x, design=NULL, tol=0.5) 
# This computes the standard deviation in the mean intensities,
# in order to compute the radius expansion required for handling
# randomly-distributed intensity shifts between samples.
#
# written by Aaron Lun
# created 27 October 2016   
{
    ci <- .raw_cellIntensities(x)
    sample.id <- .raw_sample_id(x)

    # Computing mean intensities for all (used) markers in all samples.
    nsamples <- ncol(x)
    all.means <- vector("list", nsamples)
    for (s in seq_len(nsamples)) { 
        all.means[[s]] <- rowMeans(ci[,sample.id==s,drop=FALSE])
    }
    all.means <- do.call(rbind, all.means)
    
    # Fitting a linear model to estimate variance of shift process per marker.
    if (is.null(design)) { 
        design <- matrix(1, nrow=nrow(all.means), ncol=1)
    }
    fit <- lm.fit(y=all.means, x=design)

    # Taking the mean variance across all markers, representing squared shift.
    shift2 <- mean(fit$effects[-fit$qr$pivot[seq_len(fit$qr$rank)],]^2) 
    
    # Adding to default tolerance to get per-marker radius.
    # This is done based on squared distance, hence the squaring and rooting.
    # Doubling the shift to account for distances between cells, not just distance from cell to mean.
    return(sqrt(2*shift2 + tol^2))
}
