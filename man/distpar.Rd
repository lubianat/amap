\name{distpar}
\title{Parallelized Distance Matrix Computation}
\usage{
distpar(x, method = "euclidean", nbproc = 2, diag = FALSE, upper = FALSE)

}
\alias{distpar}

\arguments{
  \item{x}{numeric matrix or (data frame).  Distances between the rows of
    \code{x} will be computed.}
  \item{method}{the distance measure to be used. This must be one of
    \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
    \code{"canberra"} or \code{"binary"} \code{"pearson"} or
    \code{"correlation"}.
    Any unambiguous substring can be given.}
  \item{nbproc}{Integer, Number of subprocess for parallelization}
  \item{diag}{logical value indicating whether the diagonal of the
    distance matrix should be printed by \code{print.dist}.}
  \item{upper}{logical value indicating whether the upper triangle of the
    distance matrix should be printed by \code{print.dist}.}
}
\description{
  This function computes and returns the distance matrix computed by
  using the specified distance measure to compute the distances between
  the rows of a data matrix.
}


\value{
  An object of class \code{"dist"}.


}

\seealso{
  \code{\link{dist}}
}
\examples{
x <- matrix(rnorm(100), nrow=5)
dist(x)
dist(x, diag = TRUE)
dist(x, upper = TRUE)
m <- as.matrix(dist(x))
d <- as.dist(m)
stopifnot(d == dist(x))
names(d) <- LETTERS[1:5]
print(d, digits = 3)

## example of binary and canberra distances.
x <- c(0, 0, 1, 1, 1, 1)
y <- c(1, 0, 1, 1, 0, 1)
dist(rbind(x,y), method="binary")
## answer 0.4 = 2/5
dist(rbind(x,y), method="canberra")
## answer 2 * (6/5)
}
\keyword{multivariate}
\keyword{cluster}
