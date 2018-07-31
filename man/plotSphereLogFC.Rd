\name{plotSphereLogFC}
\alias{plotSphereLogFC}
\alias{plotCellLogFC}

\title{Plot changes in hypersphere abundance}
\description{Visualize hyperspheres in low-dimensional space, coloured by log-fold change in abundance for each hypersphere.}

\usage{
plotSphereLogFC(x, y, logFC, max.logFC=NULL, zero.col=0.8, 
    length.out=100, pch=16, ...)

# Deprecated - all arguments passed to plotSphereLogFC.
plotCellLogFC(...)
}

\arguments{
\item{x, y}{A numeric vector of coordinates for each hypersphere.}
\item{logFC}{A numeric vector of log-fold changes for each hypersphere.}
\item{max.logFC}{A numeric scalar specifying the maximum absolute log-fold change.}
\item{zero.col}{A numeric scalar between 0 and 1, specifying the greyscale intensity to represent a log-fold change of zero.}
\item{length.out}{An integer scalar specifying the resolution of the colour bar.}
\item{pch, ...}{Additional arguments to pass to \code{plot}.}
}

\details{
Each hypersphere is represented by a point in the two-dimensional embedding, coloured from blue (negative log-FC) to grey (zero log-FC) to red (positive log-FC).
The darkness of the grey colour is set with \code{zero.col}.

If \code{max.logFC} is not \code{NULL}, extreme values in \code{logFC} are winsorized to lie within \code{[-max.logFC, max.logFC]}.
This preserves the resolution of colours for smaller log-fold changes.
}

\value{
A plot of the low-dimensional embedding of the hypersphere locations is made on the current graphics device.
A vector of colours of length \code{length.out} is returned, containing the colour gradient used for the log-fold changes. 
The vector names contains the numeric values associated with each colour.
This can be used to construct a colour bar with \code{\link{createColorBar}}.
}

\author{
Aaron Lun
}

\seealso{
\code{\link{viridis}},
\code{\link{intensityRanges}}
}

\examples{
# Making up some coordinates.
x <- rnorm(100)
y <- rnorm(100) 

# Log-FC plot and colour bar.
logFC <- rnorm(100)
out <- plotSphereLogFC(x, y, logFC)

out <- plotSphereLogFC(x, y, logFC, max.logFC=0.5)
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
createColorBar(out) 
text(-0.6, 0, srt=90, "Log-FC", cex=1.2)
}

