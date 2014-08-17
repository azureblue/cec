\name{mouseset}
\alias{mouseset}

\title{
Mouse set
}
\description{
Creates matrix of \code{dim}-dimensional points that form a "mouse-like" set with uniform density.
}
\usage{
mouseset(n = 4000, r.head = 2, r.left.ear = 1.1, r.right.ear = 1.1, left.ear.dist = 2.5,
right.ear.dist = 2.5, dim = 2)
}
\arguments{
  \item{n}{
Number of points to generate.
}
  \item{r.head}{
Radius of mouse head.
}
  \item{r.left.ear}{
Radius of mouse left ear.
}
  \item{r.right.ear}{
Radius of mouse right ear.
}
  \item{left.ear.dist}{
Distance between the center of head and the of left ear.
}
  \item{right.ear.dist}{
Distance between the center of head and the of right ear.
}
  \item{dim}{
Dimension of points.
}
}

\value{
Matrix of points with \code{n} rows and \code{dim} cols.
}

\seealso{
\code{\link{ball}}
}
\examples{
M = mouseset(n=7000, r.head=2, r.left.ear=1.1, r.right.ear=1.1, left.ear.dist=2.5, 
right.ear.dist=2.5, dim=2)
plot(M, cex = 0.5, pch = 16)
}
\keyword{ ~mouse}
\keyword{ ~mouseset}
\keyword{ ~set}
\keyword{ ~uniform}
