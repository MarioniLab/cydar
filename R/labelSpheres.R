#' @export
#' @importFrom kmknn queryKNN
labelSpheres <- function(x, labels)
# Spreads the labels to all hyperspheres, based on the closest labelled hypersphere 
# 
# written by Aaron Lun 
# created 15 May 2017
{
    if (is(x, "CyData")) {
        .check_cell_data(x)
        coords <- .raw_intensities(x)
    } else {
        coords <- x
    }

    stopifnot(identical(nrow(coords), length(labels)))
    has.label <- labels!=""
    if (!any(has.label) || all(has.label)) { 
        return(labels) 
    }

    fresh.labels <- labels[has.label]
    ulabels <- unique(fresh.labels)
    if (length(ulabels)==1) { 
        labels[] <- ulabels
        return(labels)
    }

    # Finding closest labelled hypersphere for each other hypersphere.
    coords <- as.matrix(coords)    
    labelled <- coords[has.label,,drop=FALSE]
    closest <- queryKNN(X=labelled, query=coords, k=1, get.distance=FALSE)
    fresh.labels[closest$index]
}
