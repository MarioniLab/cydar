\name{cydarnews}
\title{cydar News}
\encoding{UTF-8}

\section{Version 1.30.0}{\itemize{
\item Updated to use the latest \pkg{BiocNeighbors}, which means that the output of \code{prepareCellData()} is no longer serializable.
This is because the \code{precomputed} index is now an external pointer to a C++-owned data structure, and cannot be moved between sessions.
(Arguably, it was a mistake to expose these internals in the expected workflow for this package, but that ship has sailed.)
Hopefully, this should not affect most users as they should be primarily interacting with the \pkg{CyData} object.
}}

\section{Version 1.12.0}{\itemize{
\item \code{prepareCellData()} now returns a List rather than a half-complete CyData object.
\item \code{expandRadius()} and \code{neighborDistances()} work explicitly on the List produced by \code{prepareCellData()}.
\item \code{countCells()} will compute weighted-median intensities for unused markers as well.
Minor bugfix to avoid stochastic results due to numerical imprecision when computing medians.
\item \code{intensities()} has a new \code{mode=} argument to return intensities for unused markers.
}}

\section{Version 1.10.0}{\itemize{
\item Removed deprecated plotCell* functions.    
}}

\section{Version 1.8.0}{\itemize{
\item Modified the plotSphere*() functions for more flexibility in colour choice.

\item Added fix.zero= option in normalizeBatch() for range-based normalization.
}}

\section{Version 1.6.0}{\itemize{
\item Restructured the CyData class for simplicity and internal fields.

\item Deprecated plotCell* functions, renamed them to plotSphere*.

\item Added the createColorBar() convenience function.

\item Removed the diffIntDist() function.

\item Restored option for quantile normalization in normalizeBatch().
Switched to deterministic algorithm for sampling when mode="warp".
}}

\section{Version 1.4.0}{\itemize{
\item Bug fix to interpreSpheres() when making additional plots.

\item Switched to custom colour calculation in plotCellIntensity().
}}

\section{Version 1.2.0}{\itemize{
\item Added labelSpheres() function for labelling unannotated hyperspheres.

\item Exported multiIntHist() for plotting multiple intensity histograms.

\item Slight fix to spatialFDR(), which now computes the correct n-th nearest neighbour.
}
}

\section{Version 1.0.0}{\itemize{
\item 
New package \code{cydar}, for detecting differential abundance in mass cytometry data.
}}
