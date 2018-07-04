\name{CyData-class}
\docType{class}

\alias{CyData}
\alias{CyData-class}
\alias{show,CyData-method}

\alias{[,CyData,ANY-method}
\alias{[,CyData,ANY,ANY-method}
\alias{[,CyData,ANY,ANY,ANY-method}
\alias{[<-,CyData,ANY,ANY,CyData-method}
\alias{cbind,CyData-method}

\alias{markernames}
\alias{markernames,CyData-method}
\alias{cellAssignments}
\alias{cellAssignments,CyData-method}
\alias{intensities}
\alias{intensities,CyData-method}

\title{CyData class and methods}
\description{An overview of the CyData class and applicable methods.}

\usage{
\S4method{cellAssignments}{CyData}(x)

\S4method{intensities}{CyData}(x)

\S4method{markernames}{CyData}(object, all=FALSE)

CyData(...)
}

\arguments{
\item{x, object}{A CyData object.}
\item{all}{A logical scalar indicating whether all marker names should be returned.
By default, only the markers selected for use in \code{\link{prepareCellData}} are returned.}
\item{...}{Arguments to be passed to the \code{\link{SummarizedExperiment}} constructor.}
}

\details{
The CyData class is derived from the \linkS4class{SummarizedExperiment} class.
It is intended to store the cell counts for each group of cells (e.g., hyperspheres, clusters, in the rows) for each sample (columns).
It also stores the median intensities for each group and the identity of cells in the groups, parallel to the rows.

CyData objects should not be created directly by users, i.e., the \code{CyData} constructor should not be called manually.
It is safer to rely on functions like \code{\link{prepareCellData}} and \code{\link{countCells}} to create the objects.

The subsetting and combining behaviour of CyData objects is mostly the same as that of \linkS4class{SummarizedExperiment} objects.
The only difference is that any attempt to subset or combine by column will result in a warning.
This is because the columns are usually not independent in contexts involving clustering cells across multiple samples.
If a sample is to be removed, it is more appropriate to do so in the function that generates the CyData object (usually \code{\link{prepareCellData}}).
}

\value{
The \code{CyData} constructor returns a CyData object containing the specified information.

The \code{intensities} method returns a matrix of intensities for each group of cells (rows) and markers (columns).

The \code{cellAssignments} method returns a list of integer vectors, where each vector contains the identities of the cells in a group.
}

\examples{
# A minimal example.
my.markers <- DataFrame(row.names=LETTERS)
cyd.minimal <- CyData(markerData=my.markers)

# Adding extra detail.
counts <- matrix(rpois(1000, 10), ncol=10)
medians <- matrix(rgamma(1000, 1, 1), ncol=10)
cell.int <- matrix(rgamma(10000, 1, 1), nrow=10)
marker.data <- DataFrame(row.names=LETTERS[1:10])
cell.data <- DataFrame(sample.id=sample(10, 1000, replace=TRUE))
cyd <- CyData(assay=counts, markerData=marker.data, cellData=cell.data,
    intensities=medians, cellIntensities=cell.int)
}

\author{
Aaron Lun
}
