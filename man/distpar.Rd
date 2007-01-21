\name{distpar}
\title{Parallelized Distance Matrix Computation}
\usage{
distpar(x, method = "euclidean", nbproc = 2, diag = FALSE, upper = FALSE)

}
\alias{distpar}

\arguments{
  \item{x}{numeric matrix or (data frame) or an object of class
    "exprSet".  Distances between the rows of
    \code{x} will be computed.}
  \item{method}{the distance measure to be used. This must be one of
    \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
    \code{"canberra"}, \code{"binary"}, \code{"pearson"},
    \code{"correlation"} or \code{"spearman"}.
    Any unambiguous substring can be given.}
  \item{nbproc}{integer, Number of subprocess for parallelization}
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

\references{
  Antoine Lucas and Sylvain Jasson, \emph{Using amap and ctc Packages
    for Huge Clustering}, R News, 2006, vol 6, issue 5 pages 58-60.
  }


\seealso{
  \code{\link{Dist}}
}
\examples{
x <- matrix(rnorm(100), nrow=5)

## compute dist with 8 threads
distpar(x,nbproc=8)

## compute pearson dist with 8 threads
distpar(x,nbproc=8,method="pearson")


}
\keyword{multivariate}
\keyword{cluster}
