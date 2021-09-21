\name{plot.cec}
\alias{plot.cec}

\title{
Plot CEC.
}

\description{
Presents the results of \code{cec} function in the form of a plot. 
Colors of data points depend of the cluster they belong to. 
Ellipses are drawn with regards to the covariance (either model or sample) of each cluster.
}

\usage{
\method{plot}{cec}(x, col, cex = 0.5, pch = 19, cex.centers = 1, pch.centers = 8, 
ellipses.lwd = 4, ellipses = TRUE, model = TRUE, xlab, ylab, ...)
}

\arguments{

\item{x}{
The result of \code{cec} function.
}

\item{col}{
Use this argument to change default colors of points in the clusters.
}

\item{cex}{
Basically the size of the points, see \link{points}.
}

\item{pch}{
See \link{points}.
}

\item{cex.centers}{
The same as \code{cex} parameter, except that it's related to the centers' means.
}

\item{pch.centers}{
The same as \code{pch} parameter, except that it's related to the centers' means.
}

\item{ellipses.lwd}{
Width of ellipses, \link{points}.
}

\item{ellipses}{
If this parameter is TRUE, ellipses will be drawn.
}

\item{model}{
If this parameter is TRUE, the model (expected) covariance will be used for each cluster insted of sample covariance (MLE) of the points in the cluster, when drawing ellipses.
}

\item{xlab}{
See \link{plot}.
}

\item{ylab}{
See \link{plot}.
}

\item{\dots}{
Arguments are passed to \code{plot} function when drawing data points.
}
}

\value{No return value.}

\seealso{
\code{\link{print.cec}}
}

\examples{
## See the examples of function cec.
}

\keyword{hplot}

