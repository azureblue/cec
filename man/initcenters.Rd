\name{initcenters}
\alias{initcenters}
\title{
Center initialization
}
\description{
Creates a matrix of \code{k} points (centers) based on a given matrix of points. 
One of two method can be used: Kmeans++ centers initialization method or a random choice of data points.
}
\usage{
initcenters(x, k, method = c("kmeans++", "random"))
}
\arguments{

\item{x}{
Dataset as a matrix of n-dimensional points.
}

\item{k}{
Number of points (centers) to generate.
}

\item{method}{
Generation method. Possible values are: "kmeans++", "random.points".
}

}

\value{
Matrix points (centers) with \code{k} rows.
}

\examples{
m = matrix(runif(3000), 1000, 3)
initcenters(m, 3, method = "km")
}

\keyword{ ~centers }
\keyword{ ~initialization }
