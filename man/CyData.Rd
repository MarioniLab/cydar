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

\alias{cellIntensities}
\alias{cellInformation}
\alias{getCenterCell}

\title{CyData class and methods}
\description{An overview of the CyData class and applicable methods.}

\usage{
\S4method{cellAssignments}{CyData}(x)

\S4method{intensities}{CyData}(x)

\S4method{markernames}{CyData}(object, mode=c("used", "all", "unused"))

cellIntensities(x, mode=c("used", "all", "unused"))

cellInformation(x)

getCenterCell(x)

CyData(...)
}

\arguments{
\item{x, object}{A CyData object.}
\item{mode}{String specifying the types of markers that should be returned.
By default, only the markers used for distance calculations are returned (i.e., set as \code{markers=} in \code{?\link{prepareCellData}}).}
\item{...}{Arguments to be passed to the \code{\link{SummarizedExperiment}} constructor.}
}

\details{
The CyData class is derived from the \linkS4class{SummarizedExperiment} class.
It is intended to store the cell counts for each group of cells (rows) for each sample (columns).
Groups are intended to be hyperspheres (see \code{\link{countCells}}) but could also be arbitrary clusters of cells.
It also stores the median intensities for each group and the identity of cells in the groups, parallel to the rows.

CyData objects should not be created directly by users, i.e., the \code{CyData} constructor should not be called manually.
The class has some strict validity conditions that are not easily satisfied by manual construction.
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

The \code{markernames} method returns a character vector of marker names, by default corresponding to the columns of the \code{intensities}.
If \code{mode="unused"}, names of all unused markers are returned (i.e., not set as \code{markers=} in \code{?\link{prepareCellData}}).
If \code{mode="all"}, names of all markers are returned.

The \code{cellIntensities} function returns a matrix of cell-level intensities for the specified markers.
Rows correspond to markers while columns correspond to individual cells (not groups of cells, unlike \code{intensities}).

The \code{cellInformation} function returns a DataFrame of cell-level information, where each row corresponds to a column of the matrix from \code{cellIntensities}.
It contains \code{sample}, an integer field specifying the sample of origin for each cell; 
and \code{row}, an integer field specifying the row of the sample-specific intensity matrix for that cell.

The \code{getCenterCell} function returns the index of the cell used at the center of each hypersphere.
These indices correspond to columns of the \code{cellIntensities} output and rows of the \code{cellInformation} output.
}

\examples{
example(countCells, echo=FALSE)

markernames(cnt)
head(intensities(cnt))
head(cellAssignments(cnt))
}

\author{
Aaron Lun
}
