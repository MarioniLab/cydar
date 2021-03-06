\name{plotSphereIntensity}
\alias{plotSphereIntensity}

\title{Plot cell or hypersphere data}
\description{Visualize cells or hyperspheres in low-dimensional space, coloured by marker intensities or log-fold changes.}

\usage{
plotSphereIntensity(x, y, intensity, irange=NULL, 
    col.range=viridis(100), pch=16, ...)
}

\arguments{
\item{x, y}{A numeric vector of coordinates for each hypersphere.}
\item{intensity}{A numeric vector specifying the marker intensities for each hypersphere.}
\item{irange}{A numeric vector of length 2, specifying the upper and lower bound for the intensities.}
\item{col.range}{A vector of colours specifying the colour scale to be used for increasing intensity.
More values represent a higher-resolution scale.}
\item{pch, ...}{Additional arguments to pass to \code{plot}.}
}

\details{
Each hypersphere is represented by a point in the two-dimensional embedding, colored using the viridis colour scheme, i.e., purple (low intensity) to green (medium) to yellow (high).
If \code{irange} is not \code{NULL}, extreme values in \code{intensity} will be winsorized to lie within \code{irange}.
This preserves the resolution of colours for smaller changes at low intensities. 
Users should consider using \code{\link{intensityRanges}} to define appropriate values of \code{irange} for each marker.
}

\value{
A plot of the low-dimensional embedding of the hypersphere locations is made on the current graphics device.
A vector of colours equal to \code{col.range} is returned, containing the colour gradient used for the intensities.
The vector names contain the numeric values associated with each colour.
This can be used to construct a colour bar with \code{\link{createColorBar}}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{viridis}},
\code{\link{intensityRanges}},
\code{\link{createColorBar}}
}

\examples{
# Making up some coordinates.
x <- rnorm(100)
y <- rnorm(100) 

# Intensity plot and colour bar.
intensities <- rgamma(100, 2, 2)
out <- plotSphereIntensity(x, y, intensities)

plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
createColorBar(out) 
text(-0.6, 0, srt=90, "Intensity", cex=1.2)
}

