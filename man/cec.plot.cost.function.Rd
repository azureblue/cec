\name{cec.plot.cost.function}
\alias{cec.plot.cost.function}

\title{
Plot cost function of CEC results.
}

\description{
Presents cost function changes at each iteration of the CEC algorithm in the form of plot.
}

\usage{
cec.plot.cost.function(C, xlab = "Iteration", ylab = "Cost function", lwd = 5, 
col = "red", lwd.points = 5, pch.points = 19, col.points = "black", 
mgp = c(1.5, 0.5, 0), ...)
}

\arguments{

\item{C}{
The result of \code{cec} function.
}

\item{xlab}{
See \link{title}.
}

\item{ylab}{
See \link{title}.
}

\item{lwd}{
Width of lines. See \link{points}.
}

\item{col}{
Color of lines. See \link{points}.
}

\item{lwd.points}{
Size of points. See \link{points}.
}

\item{pch.points}{
Plotting 'character' of points. See \link{points}.
}

\item{col.points}{
Color of points.See \link{points}.
}

\item{mgp}{
See \link{par}.
}

\item{\dots}{
Passed to line plot function.
}

}
\seealso{
\link{plot.cec}
}
\examples{
## See the examples of function cec.
}

\keyword{ ~plot }
\keyword{ ~cec }
