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
Either a matrix of initial centers or a the number of initial centers (\code{k}). 
In the latter case, initial centers will be generated using a method depending on
the \code{centers.init} argument.
}

\item{iter.max}{
Maximum number of iterations at each start.
}

\item{nstart}{
Number of clusterings to perform (with different initial centers). Only the best 
clustering (with the lowest energy) will be returned. Value grater then one is valid 
only if the \code{centers} argument is a number.
}

\item{centers.init}{
Centers initialization method. Possible values are: "kmeans++", "random".
}

\item{type}{
Type (or types) of clustering (density family). This can be either a single value 
or a vector of length equal to the number of centers. Possible values are: "covariance",
"fixedr", "spherical", "diagonal", "all".
}

\item{param}{
Parameter (or parameters) specific to a particular type of clustering. 
Not all types of clustering require parameter. Types that require parameter: 
"covariance" (matrix parameter), "fixedr" (numeric parameter). 
This can be a vector or a list (when one of the parameters is a matrix). 
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
Interactive mode. If TRUE, result of clustering will be plotted after every iteration.
}

\item{readline}{
Used only in the interactive mode. If \code{readline} is TRUE, at each iteration, before
plotting it will wait for the user to press <Return> instead of standard "before plotting"
(\code{par(ask = TRUE)}) waiting.
}

}
\value{
  Returns an object of class "cec" with available components:
  "data", "cluster", "centers", "iterations", "nclusters", "energy", "covariances", 
  "covariances.model", "time".
}

\details{
General idea of Cross Entropy Clustering relies on finding the splitting of \eqn{X \subset \R^N}{X⊂Rᴺ} into pairwise disjoint sets \eqn{X_1,\ldots,X_k}{X1, \ldots, Xk} such that the overall energy of the clusters is minimal.

The CEC algorithm minimizes the cost function which is a small modification of the MLE (Maximum Likelihood Estimation): \deqn{CEC(f,X):= -\sum_{i=1}^m \ln(\max(p_1 f_1(x_i),\ldots,p_k f_k(x_i)))}{CEC(f,X):= -\sumln(max(p1 f1(xi), \ldots, pk fk(xi)))}
where \eqn{p_i \in [0,1], \sum p_i=1}{0 ≤ pi ≤ 1, \sumpi=1  } and \eqn{f_1, \ldots, f_k}{f1, ..., fk} are Gaussian densities.

Five density models (types of clustering) have been proposed: 
\itemize{
\item"covariance" - Gaussian densities with a given covariance. The shapes of clusters depend on the given covariance 
matrix (additional parameter).
\item"fixedr" - Special case of "covariance", where the covariance matrix equals \emph{rI} for the given \emph{r} 
(additional parameter). The clustering will have a tendency to divide data into balls with approximate radius of the square root of \emph{r}.
\item"spherical" - Spherical Gaussian densities (covariance proportional to the indentity). Clusters will have a tendency 
to form balls of arbitrary sizes.
\item"diagonal" - Gaussian densities with diagonal covariane. Data will form ellipsoids with radiuses parallel to coordinate axes.
\item"all" - All Gaussian densities. Data will form ellipsoids with arbitrary radiuses. 
}

The implementation of \code{cec} function allows mixing clustering types.
}
\author{
Konrad Kamieniecki, Jacek Tabor, Przemysław Spurek
}

\seealso{
\code{\link{cec-package}.}
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