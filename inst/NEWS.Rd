\name{cydarnews}
\title{cydar News}
\encoding{UTF-8}

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
