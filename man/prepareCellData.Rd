\name{prepareCellData}
\alias{prepareCellData}

\title{Prepare mass cytometry data}

\description{Convert single-cell marker intensities from a mass cytometry experiment into a format for efficient counting.}

\usage{
prepareCellData(x, markers=NULL, ...)
}

\arguments{
\item{x}{A named list of numeric matrices, where each matrix corresponds to a sample and contains expression intensities for each cell (row) and each marker (column).
Alternatively, a ncdfFlowSet object containing the same information.}
\item{markers}{A vector specifying the markers to use in downstream analyses.}
\item{...}{Additional arguments to pass to \code{\link{buildKmknn}}.}
}

\details{
This function constructs a CyData object from the marker intensities of each cell in one or more samples.
The \code{\link{buildKmknn}} function is used to precompute internal structures for downstream nearest-neighbour searching.

If \code{markers} is specified, only the selected markers will be used in the precomputation.
This restricts the markers that are used in downstream functions - namely, \code{\link{countCells}} and \code{\link{neighborDistances}}.
By default, \code{markers=NULL} which means that all supplied markers will be used.

Markers that are \emph{not} in \code{markers} will be ignored in distance calculations.
However, their intensities are still stored in the output object, for use in functions like \code{\link{medIntensities}}.
}

\value{
A CyData object with no rows and the number of columns equal to the number of samples in \code{x}.
Sample names are stored as the column names of the output object.
Precomputed values for internal use are stored in the \code{"cydar"} field of the \code{metadata} slot.
}

\author{
Aaron Lun
}

\examples{
### Mocking up some data: ###
nmarkers <- 20
marker.names <- paste0("X", seq_len(nmarkers))
nsamples <- 8
sample.names <- paste0("Y", seq_len(nsamples))

x <- list()
for (i in sample.names) {
    ex <- matrix(rgamma(nmarkers*1000, 2, 2), ncol=nmarkers, nrow=1000)
    colnames(ex) <- marker.names
    x[[i]] <- ex
}

### Running the function: ###
cd <- prepareCellData(x)
cd
}

\seealso{
\code{\link{countCells}},
\code{\link{neighborDistances}}
}
