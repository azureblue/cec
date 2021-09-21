\name{ball}
\alias{ball}

\title{
Ball
}

\description{
Generates points that form a ball with uniform density.
}

\usage{
ball(n, r, dim)
}

\arguments{

\item{n}{
Number of points to generate.
}

\item{r}{
Radius of the ball.
}

\item{dim}{
Dimension of the points.
}

}
\value{
Matrix of points with \code{n} rows and \code{dim} cols.
}

\seealso{
\code{\link{mouseset}}
}

\examples{
M = ball(4000, 0.9)
plot(M, cex = 0.5, pch = 19)
}

\keyword{datagen}
