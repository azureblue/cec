\name{cec}
\alias{cec}
\title{
Cross-Entropy Clustering
}
\description{
Performs Cross-Entropy Clustering on a data matrix.
}
\usage{
cec(x, centers, iter.max = 20, nstart = 1, centers.init  = c("kmeans++", "random"), 
type = c("covariance", "fixedr", "spherical", "diagonal", "all"), param, 
card.min = "5\%", keep.removed = F, interactive = F, readline = T)
}
\arguments{
  \item{x}{
Numeric matrix of data. 
}
  \item{centers}{
Either a matrix of initial centers or the number of initial centers (\code{k}). 
In the latter case, initial centers will be generated using method depending on
\code{centers.init} argument.
}
  \item{iter.max}{
Maximum number of iterations at each start.
}
  \item{nstart}{
Number of clusterings to perform (with different iniial centers). Only the best 
clustering (with the lowest energy) will be returned. Value grater then one is valid 
only if \code{centers} argument is a number.
}
  \item{centers.init}{
Centers initialization method. Possible values are: "kmeans++", "random".
}
  \item{type}{
Type (or types) of the clustering (density family). This can be either a single value 
or a vector of length equal to the number of centers. Possible values are: "covariance",
"fixedr", "spherical", "diagonal", "all".
}

\item{param}{
Parameter (or parameters) specific to particular type of the clustering. Not all types of
clustering require parameter. This can be a vector or a list (when one of the parameter is a matrix). 
}

\item{card.min}{
Minimal cluster cardinality. If cluster cardinality becomes less than card.min, cluster
is removed. This argument can be either integer number or string ended with percent sign
(e.g. "5\%").
}

\item{keep.removed}{
If this parameter is FASLE, removed clusters will be visible in the results as NA in centers matrix (as well as in the list of covariances).
}

\item{interactive}{
Interactive mode. Plotting after every iteration if this arg is TRUE.
}
\item{readline}{
Used only in interactive mode. If \code{readline} is TRUE, at each iteration, before
plotting it waits for the user to press <Return> instead of standard "before plotting"
(\code{par(ask = TRUE)}) waiting.
}

}
\value{
  Returns object of class "cec" with available components:
  "data", "cluster", "centers", "iterations", "clusters.number", "energy", "covariances",
  "time".
}
\author{
Konrad Kamieniecki
}
\references{
Spurek, P. and Tabor, J. (2012)
Cross-Entropy Clustering
\emph{arXiv preprint} arXiv:1210.5594
}

\examples{
# Cross-Entropy Clustering
#
#
# 
#
cat("Example of clustering random dataset of 3 Gaussians using spherical Gaussian densities,
10 random initial centers and 7\% as minimal cluster size.")
m1 = matrix(rnorm(2000, sd=1), ncol=2)
m2 = matrix(rnorm(2000, mean = 3, sd = 1.5), ncol = 2)
m3 = matrix(rnorm(2000, mean = 3, sd = 1), ncol = 2)
m3[,2] = m3[,2] - 5
m = rbind(m1, m2, m3)
centers = initcenters(m, 10)
par(ask = TRUE)
plot(m, cex = 0.5, pch = 16)
Z = cec(m, centers, type="sp", iter.max = -1, card.min="7\%")
plot(Z)
Z = cec(m, centers, type="sp", iter.max = 100, card.min="7\%")
plot(Z)
# Finished
#
# Result:
Z
cec.plot.energy(Z)
cat ("Example of clustering mouse-like set. Press <Enter>.")
m = mouseset(n=7000, r.head=2, r.left.ear=1.1, r.right.ear=1.1, left.ear.dist=2.5,
right.ear.dist=2.5, dim=2)
plot(m, cex = 0.5, pch = 16)
centers = initcenters(m, 3)
Z = cec(m, centers, type="sp", iter.max = -1, card.min="5\%")
plot(Z)
Z = cec(m, centers, type="sp", iter.max = 100, nstart=3, card.min="5\%")
plot(Z)
# Finished
#
# Result:
Z
cec.plot.energy(Z)
cat ("Example of clustering uniform density ball starting with 2 clusters")
m = ball(7000, 1)
plot(m, cex = 0.5, pch = 16)
centers = initcenters(m, 2)
Z = cec(m, centers, type="sp", iter.max = -1, card.min="5\%")
plot(Z)
Z = cec(m, centers, type="sp", iter.max = 100, card.min="5\%")
plot(Z)
# Finished
#
# Result:
Z
cec.plot.energy(Z)
}
\keyword{ ~cec }
\keyword{ ~clustering }
