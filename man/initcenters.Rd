\name{initcenters}
\alias{initcenters}
\title{
Center initialization
}
\description{
Creates matrix of points (centers) based on given matrix of points. One of two method may be used: Kmeans++ centers initialization method or random choice of data points.
}
\usage{
initcenters(x, c, method = c("kmeans++", "random"))
}
\arguments{
  \item{x}{
Dataset as a matrix of n-dimensional points.
}
  \item{c}{
Number of centers (points) to be computed.
}
  \item{method}{
Method of computing centers. Possible values are: "kmeans++", "random.points".
}
}
\value{
N-cols matrix representing set of points (centers).
}

\examples{
m = matrix(runif(3000), 1000, 3)
initcenters(m, 3, method = "km")
}

\keyword{ ~centers }
\keyword{ ~initialization }
