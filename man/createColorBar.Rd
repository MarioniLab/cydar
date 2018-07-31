\name{createColorBar}
\alias{createColorBar}

\title{Create a color bar}
\description{A convenience function to create a color bar on the current graphics device.}

\usage{
createColorBar(colors, top.text=NULL, bottom.text=NULL, 
    lower=-0.5, upper=0.5, x.pos=0, width=1, cex=1.5) 
}

\arguments{
    \item{colors}{A character vector of colors to be used in successive steps of the color bar, from bottom to top.}
    \item{top.text}{String containing the label to be placed at the top of the color bar.}
    \item{bottom.text}{String containing the label to be placed at the bottom of the color bar.}
    \item{lower}{Numeric scalar specifying the bottom y-coordinate of the color bar.}
    \item{upper}{Numeric scalar specifying the top y-coordinate of the color bar.}
    \item{x.pos}{Numeric scalar specifying the x-coordinate of the middle of the color bar.}
    \item{width}{Numeric scalar specifying the width of the color bar.}
    \item{cex}{Numeric scalar specifying the size of the text.}
}

\details{
If \code{bottom.text} or \code{top.text} are not specified, they will be take from the first and last names of \code{colors}, respectively.
If no names are available, no labels will be created.
}

\value{
A color bar is created on the current graphics device.
\code{NULL} is invisibly returned.
}

\author{
Aaron Lun
}

\examples{
plot(0,0, axes=FALSE, xlab="", ylab="")
createColorBar(heat.colors(100), "High", "Low")
}
