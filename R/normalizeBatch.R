#' @export
normalizeBatch <- function(batch.x, batch.comp, mode="range", p=0.01, target=NULL, markers=NULL, ...)
# Performs warp- or range-based adjustment of different batches, given a 
# list of 'x' objects like that used for 'prepareCellData'
# and another list specifying the composition of samples per batch.
#
# written by Aaron Lun
# created 27 October 2016
{
    if (is.null(batch.comp)) {
        batch.comp <- lapply(batch.x, function(i) rep(1, length(i)))
    }
    nbatches <- length(batch.x)
    if (nbatches!=length(batch.comp)) {
        stop("length of 'batch.x' and 'batch.comp' must be identical")
    }

    # Checking the number of markers we're dealing with.
    batch.out <- vector("list", nbatches)
    for (b in seq_len(nbatches)) { 
        out <- .pull_out_data(batch.x[[b]])

        if (!is.null(markers)) {
            mm <- match(markers, out$markers)
            if (any(is.na(mm))) { stop("some 'markers' not present in batch") }
            out$markers <- out$markers[mm]
            out$exprs <- lapply(out$exprs, function(x) { x[,mm,drop=FALSE] })
        }
        batch.out[[b]] <- out

        if (b==1L) {
            ref.markers <- out$markers
        } else if (!identical(ref.markers, out$markers)) { 
            stop("markers are not identical between batches")
        }
        if (length(out$samples)!=length(batch.comp[[b]])) {
            stop("corresponding elements of 'batch.comp' and 'batch.x' must have same lengths")
        }
    }

    # Expanding possible modes.    
    mode <- rep(mode, length.out=length(ref.markers))
    if (is.null(names(mode))) {
        names(mode) <- ref.markers
    }

    # Checking 'target' specification.
    if (!is.null(target)) { 
        target <- as.integer(target)
        if (target < 1L || target > nbatches) {
            stop("'target' must be a positive integer no greater than the number of batches")
        }
    }

    # Calculating weights.
    batch.weights <- .computeCellWeights(batch.out, batch.comp)
    EXTRACTOR <- .pullOutMarkers(batch.out, batch.weights)

    # Setting up an output object.
    output <- vector("list", nbatches)
    for (b in seq_len(nbatches)) { 
        cur.out <- batch.out[[b]]
        nsamples <- length(cur.out$samples)

        cur.exprs <- vector("list", nsamples)
        for (s in seq_len(nsamples)) {
            cur.exprs[[s]] <- cur.out$exprs[[s]]
            colnames(cur.exprs[[s]]) <- cur.out$markers
        }
        names(cur.exprs) <- cur.out$samples
        output[[b]] <- cur.exprs 
    }
    names(output) <- names(batch.x)

    for (m in ref.markers) {
        for.norm <- EXTRACTOR(m)
        all.obs <- for.norm$exprs
        all.wts <- for.norm$weights
        
        # Choosing the normalization method.
        curmode <- match.arg(mode[m], c("none", "range", "warp", "quantile"))
        if (curmode=="none") {
            ;
        } else if (curmode=="warp") { 
            converters <- .transformDistr(all.obs, all.wts, m, target=target, ...)
        } else if (curmode=="range") {
            converters <- .rescaleDistr(all.obs, all.wts, target=target, p=p)
        } else if (curmode=="quantile") {
            converters <- .quantileDistr(all.obs, all.wts, target=target)
        }

        # Applying the normalization method.
        for (b in seq_len(nbatches)) { 
            converter <- converters[[b]]
            cur.out <- batch.out[[b]]
            for (s in seq_along(cur.out$exprs)) {
                output[[b]][[s]][,m] <- converter(cur.out$exprs[[s]][,m])                
            }
        }        
    }
    return(output)
}

############################################################

.computeCellWeights <- function(batch.out, batch.comp) 
# Estimates the weight for each sample in each group in each batch.
# This accounts for differences in numbers of groups per batch,
# and for differences in numbers of cells per sample.
{ 
    all.levels <- unique(unlist(batch.comp))
    batch.comp <- lapply(batch.comp, factor, levels=all.levels)

    comp.batches <- do.call(rbind, lapply(batch.comp, table))
    ref.comp <- colMeans(comp.batches)
    batch.weight <- t(ref.comp/t(comp.batches))

    empty.factors <- colSums(!is.finite(batch.weight)) > 0
    if (all(empty.factors)) {
        stop("no level of 'batch.comp' is common to all batches")
    }
    batch.weight <- batch.weight[,!empty.factors,drop=FALSE]

    nbatches <- length(batch.out)
    batch.weights <- vector("list", nbatches)
    for (b in seq_len(nbatches)) { 
        cur.comp <- batch.comp[[b]]
        cur.out <- batch.out[[b]]
        cur.weights <- numeric(length(cur.out$exprs))
        
        for (s in seq_along(cur.out$exprs)) {
            sample.level <- as.character(cur.comp[s])
            num.cells <- nrow(cur.out$exprs[[s]])
            if (sample.level %in% colnames(batch.weight)) { 
                cur.weights[s] <- 1/num.cells * batch.weight[b,sample.level]
            } else {
                cur.weights[s] <- NA_real_
            }
        }
        batch.weights[[b]] <- cur.weights
    }
    return(batch.weights)
}

.pullOutMarkers <- function(batch.out, batch.weights) 
# Returns functions that pull out intensities and weights for every marker.
# Weights are also returned here (but not recalculated) to synchronize removal of zero-weight samples.
{
    nbatches <- length(batch.out)
    all.mats <- all.wts <- vector("list", nbatches)
    for (b in seq_along(batch.weights)) {
        cur.weights <- batch.weights[[b]]
        keep <- !is.na(cur.weights)
        all.mats[[b]] <- batch.out[[b]]$exprs[keep]
        num.cells <- vapply(all.mats[[b]], nrow, FUN.VALUE=0L)
        all.wts[[b]] <- rep(cur.weights[keep], num.cells)
    }

    function(m) { 
        all.obs <- vector("list", nbatches)
        for (b in seq_len(nbatches)) { 
            cur.out <- all.mats[[b]]
            nsamples <- length(cur.out)

            cur.obs <- vector("list", nsamples)
            for (s in seq_len(nsamples)) { 
                cur.obs[[s]] <- cur.out[[s]][,m]
            }
            all.obs[[b]] <- unlist(cur.obs)
        }
        list(exprs=all.obs, weights=all.wts)
    }
}

############################################################

.getECDF <- function(cur.obs, cur.wts) 
# Computes the bits and pieces necessary for the ECDF.
{
    o <- order(cur.obs)
    cur.obs <- cur.obs[o]
    cur.wts <- cur.wts[o]

    # Taking the midpoint of each step, rather than the start/end points. 
    mid.cum.weight <- cumsum(cur.wts) - cur.wts/2
    total.weight <- sum(cur.wts)
    list(cumprob=mid.cum.weight/total.weight, value=cur.obs)
}

#' @importClassesFrom flowCore flowSet
#' @importFrom flowCore flowFrame normalization
#' @importFrom stats splinefun approx
#' @importMethodsFrom BiocGenerics colnames<-
#' @importMethodsFrom flowCore normalize
.transformDistr <- function(all.obs, all.wts, name, target, ...) {
    # Create a mock sample by sampling at points along the ECDF.
    nbatch <- length(all.obs)
    cur.ffs <- vector("list", nbatch)
    for (b in seq_len(nbatch)) {
        ecdf.out <- .getECDF(all.obs[[b]], all.wts[[b]])
        pts <- seq(0, 1, length.out=(length(all.obs[[b]])+2L)) # Note, always have the first and last entries.
        mock <- approx(ecdf.out$cumprob, ecdf.out$value, xout=pts, rule=2)$y
        cur.ffs[[b]] <- flowFrame(cbind(M=mock))
    }        

    if (!is.null(target)) {
        names(cur.ffs) <- seq_len(nbatch)
        target <- names(cur.ffs)[target]
    }
    fs <- as(cur.ffs, "flowSet")
    colnames(fs) <- name

    # Applying warping normalization, as described in the flowStats vignette.
    new.fs <- flowStats::warpSet(fs, name, monwrd=TRUE, target=target, ...)

    # Defining warp functions (setting warpFuns doesn't really work, for some reason).
    converter <- vector("list", nbatch)
    for (b in seq_len(nbatch)) {
        old.i <- exprs(fs[[b]])[,1]
        new.i <- exprs(new.fs[[b]])[,1]
        converter[[b]] <- splinefun(old.i, new.i)
    }
    return(converter)
}

#' @importFrom stats lm approx
.rescaleDistr <- function(all.obs, all.wts, target, p) {
    # Computing the average max/min (robustly).
    nbatches <- length(all.obs)
    batch.min <- batch.max <- numeric(nbatches)
    for (b in seq_len(nbatches)) { 
        ecdf.out <- .getECDF(all.obs[[b]], all.wts[[b]])
        out <- approx(ecdf.out$cumprob, ecdf.out$value, xout=c(p, 1-p), rule=2)$y
        batch.min[b] <- out[1]
        batch.max[b] <- out[2]
    }

    # Selecting the target batch to perform the normalization.
    if (is.null(target)) { 
        targets <- c(mean(batch.min), mean(batch.max))
    } else {
        targets <- c(batch.min[target], batch.max[target])
    }

    # Scaling intensities per batch so that the observed range equals the average range.
    converters <- vector("list", nbatches)
    FUNGEN <- function(fit) {
        m <- coef(fit)[2]
        b <- coef(fit)[1]
        function(x) { x * m + b }
    }

    for (b in seq_len(nbatches)) {
        current <- c(batch.min[b], batch.max[b])
        fit <- lm(targets ~ current)
        converters[[b]] <- FUNGEN(fit)
    }
    return(converters)
}

#' @importFrom stats approxfun
.quantileDistr <- function(all.obs, all.wts, target) 
# Performs quantile normalization across batches.
{
    nbatches <- length(all.obs)
    all.x <- all.y <- vector("list", nbatches)
    for (b in seq_len(nbatches)) { 
        ecdf.out <- .getECDF(all.obs[[b]], all.wts[[b]])
        all.x[[b]] <- ecdf.out$cumprob
        all.y[[b]] <- ecdf.out$value
    }

    all.quanfun <- mapply(approxfun, x=all.x, y=all.y, MoreArgs=list(rule=2))
    output.fun <- vector("list", nbatches)

    for (b in seq_len(nbatches)) { 
        current.probs <- all.x[[b]]
        current.quants <- all.y[[b]]
                
        # Correcting intensities for each sample in each batch; first by
        # computing the average distribution to which all others should be squeezed.
        target.profile <- 0
        if (is.null(target)) { 
            for (b2 in seq_len(nbatches)) { 
                if (b==b2) { 
                    target.profile <- target.profile + current.quants
                } else{
                    target.profile <- target.profile + all.quanfun[[b2]](current.probs)
                }
            }
            target.profile <- target.profile/nbatches
        } else {
            target.profile <- all.quanfun[[target]](current.probs)
        }
                
        # Constructing a function to do that squeezing (we could use the ordering to 
        # get exactly which entry corresponds to which original observation, but that
        # doesn't allow for samples that weren't used, e.g., due to non-common levels).
        output.fun[[b]] <- approxfun(current.quants, target.profile, rule=2) 
    }

    return(output.fun)
}
