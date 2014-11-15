\name{cec}
\alias{cec}

\title{
Cross-Entropy Clustering
}

\description{
Performs Cross-Entropy Clustering on a data matrix.
}

\usage{
cec(x, centers, type = c("covariance", "fixedr", "spherical", "diagonal",  
"eigenvalues", "all"), iter.max = 25, nstart = 1, param, 
centers.init = c("kmeans++", "random"), card.min = "5\%", keep.removed = F, 
interactive = F, readline = T)
}

\arguments{

\item{x}{
Numeric matrix of data. 
}

\item{centers}{
Either a matrix of initial centers or the number of initial centers (\code{k}, single 
number or a vector of numbers for variable number of initial centers at each start). 
In the latter case, initial centers will be generated using a method depending on
the \code{centers.init} argument.
}

\item{type}{
Type (or types) of clustering (density family). This can be either a single value 
or a vector of length equal to the number of centers. Possible values are: "covariance",
"fixedr", "spherical", "diagonal", "eigenvalues", "all" (default).
}

\item{iter.max}{
Maximum number of iterations at each start.
}

\item{nstart}{
Number of clusterings to perform (with different initial centers). Only the best 
clustering (with the lowest cost) will be returned. Value grater then one is valid 
only if the \code{centers} argument is a number.
}

\item{centers.init}{
Centers initialization method. Possible values are: "kmeans++" (default), "random".
}

\item{param}{
Parameter (or parameters) specific to a particular type of clustering. 
Not all types of clustering require parameter. Types that require parameter: 
"covariance" (matrix parameter), "fixedr" (numeric parameter), 
"eigenvalues" (vector parameter). This can be a vector or a list 
(when one of the parameters is a matrix or a vector). 
}

\item{card.min}{
Minimal cluster cardinality. If cluster cardinality becomes less than card.min, cluster
is removed. This argument can be either an integer number or a string ended with a percent 
sign (e.g. "5\%").
}

\item{keep.removed}{
If this parameter is TRUE, removed clusters will be visible in the results as NA in centers matrix 
(as well as corresponding values in the list of covariances).
}

\item{interactive}{
Interactive mode. If TRUE, the result of clustering will be plotted after every iteration.
}

\item{readline}{
Used only in the interactive mode. If \code{readline} is TRUE, at each iteration, before
plotting it will wait for the user to press <Return> instead of standard "before plotting"
(\code{par(ask = TRUE)}) waiting.
}

}
\value{
  Returns an object of class "cec" with available components:
  "data", "cluster", "probabilities", "centers", "cost.function", "nclusters", "final.cost.function", 
  "final.nclusters", "iterations", "cost", "covariances", 
  "covariances.model", "time".
}

\details{
In the context of implementation, Cross-Entropy Clustering (CEC) aims to partition \emph{m} points into \emph{k} 
clusters so as to minimize the cost function (energy \emph{\strong{E}} of the clustering) by switching the points between clusters. 
The presented method is based on the adapted Hartigan approach, where we reduce clusters which cardinalities decreased below some small prefixed level.

The energy function \emph{\strong{E}} is given by:
\deqn{E(Y_1,\mathcal{F}_1;...;Y_k,\mathcal{F}_k) = \sum\limits_{i=1}^{k} p(Y_i) \cdot (-ln(p(Y_i)) + H^{\times}(Y_i\|\mathcal{F}_i))}{
E(Y1, F1; ...; Yk, Fk) = \sum(p(Yi) * (-ln(p(Yi)) + H(Yi | Fi)))}
where \emph{Yi} denotes the \emph{i}-th cluster, \emph{p(Yi)} is the ratio of the number of points in \emph{i}-th cluster to the total number points, \emph{\strong{H}(Yi|Fi)} is the value of cross-entropy, which represents the internal cluster energy function of data \emph{Yi} defined with respect to a certain Gaussian density family \emph{Fi}, which encodes the type of clustering we consider.

The value of the internal energy function \emph{\strong{H}} depends only on the covariance matrix of the points in the cluster (computed using maximum-likelihood method). Five implementations of \emph{\strong{H}} have been proposed (expressed as a type - model - of the clustering):

\itemize{
\item"covariance" - Gaussian densities with a fixed given covariance. The shapes of clusters depend on the given covariance 
matrix (additional parameter).
\item"fixedr" - Special case of "covariance", where the covariance matrix equals \emph{rI} for the given \emph{r} 
(additional parameter). The clustering will have a tendency to divide data into balls with approximate radius proportional to the square 
root of \emph{r}.
\item"spherical" - Spherical (radial) Gaussian densities (covariance proportional to the indentity). 
Clusters will have a tendency to form balls of arbitrary sizes.
\item"diagonal" - Gaussian densities with diagonal covariane. Data will form ellipsoids with radiuses parallel to the coordinate axes.
\item"eigenvalues" - Gaussian densities with covariance matrix having fixed eigenvalues (additional parameter). The clustering will try to divide the data into ellipsoids with fixed shape rotated by an arbitrary angle. 
\item"all" - All Gaussian densities. Data will form ellipsoids with arbitrary radiuses. 
}

The implementation of \code{cec} function allows mixing of clustering types.
}
\author{
Konrad Kamieniecki, Jacek Tabor, Przemysław Spurek
}

\seealso{
\code{\link{CEC-package}.}
}
\examples{
#
#
#
# Cross-Entropy Clustering
#
#
# 

## Example of clustering random dataset of 3 Gaussians using spherical Gaussian densities, 
## 10 random initial centers and 7\% as minimal cluster size.

m1 = matrix(rnorm(2000, sd=1), ncol=2)
m2 = matrix(rnorm(2000, mean = 3, sd = 1.5), ncol = 2)
m3 = matrix(rnorm(2000, mean = 3, sd = 1), ncol = 2)
m3[,2] = m3[,2] - 5
m = rbind(m1, m2, m3)
centers = initcenters(m, 10)
par(ask = TRUE)
plot(m, cex = 0.5, pch = 16)
## Initial centers:
Z = cec(m, centers, type="sp", iter.max = -1, card.min="7\%")
plot(Z)
## Clustering result:
Z = cec(m, centers, type="sp", iter.max = 100, card.min="7\%")
plot(Z)
# Result:
Z
# Cost function:
cec.plot.cost.function(Z)
## Example of clustering mouse-like set.
m = mouseset(n=7000, r.head=2, r.left.ear=1.1, r.right.ear=1.1, left.ear.dist=2.5,
right.ear.dist=2.5, dim=2)
plot(m, cex = 0.5, pch = 16)
centers = initcenters(m, 3)
## Initial centers:
Z = cec(m, centers, type="sp", iter.max = -1, card.min="5\%")
plot(Z)
## Clustering result:
Z = cec(m, centers, type="sp", iter.max = 100, nstart=4, card.min="5\%")
plot(Z)
# Result:
Z
# Cost function:
cec.plot.cost.function(Z)
## Example of clustering data set "Tset" using "eigenvalues" clustering type.
data(Tset)
plot(Tset, cex = 0.5, pch = 16)
centers = initcenters(Tset, 2)
## Initial centers:
Z <- cec(Tset, 5, type="eigenvalues", param=c(0.02,0.002), iter.max= -1)
plot(Z)
## Clustering result:
Z <- cec(Tset, 5, "eigenvalues", param=c(0.02,0.002), nstart=4)
plot(Z)
# Result:
Z
# Cost function:
cec.plot.cost.function(Z)
}
\references{
Spurek, P. and Tabor, J. (2014)
Cross-Entropy Clustering
\emph{Pattern Recognition} \bold{47, 9} 3046--3059
}
\keyword{ ~cec }
\keyword{ ~clustering }
