\name{hc2Newick}
\alias{hc2Newick}
\title{Convert hclust objects to Newick format files}
\description{
  Convert hclust objects to Newick format files.
}
\usage{
  hc2Newick(hc, flat=TRUE)
}
\arguments{
  \item{hc}{a \code{hclust} object (as returned by the function \code{hclust}
    in the package \code{stats})}
  \item{flat}{a boolean (see section value).}
}
\value{
  If \code{flat=TRUE} the result is a string (that you can write in a
  file).

  If \code{flat=FALSE} the result is a list (of lists). Each list
  is consituted of the elements \code{left}, \code{right} and \code{dist}.  
}

\examples{
data(USArrests)
h = hclust(dist(USArrests))
write(hc2Newick(h),file='hclust.newick')
}


\author{
  Laurent (laurent@cbs.dtu.dk)
}
\keyword{manip}
\keyword{cluster}
