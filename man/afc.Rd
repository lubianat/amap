\name{afc}
\alias{afc}
\title{Correspondance factorial analysis.}

\usage{
afc(x)
}

\arguments{
  \item{x}{A contingency table, or a result of function \code{burt} or
    \code{matlogic}}
}

\description{
  Compute an acp on a contingency table tacking into account weight
  of rows and columns
}

\examples{
\dontrun{
color <- as.factor(c('blue','red','red','blue','red'))
size <- as.factor(c('large','large','small','medium','large'))
x <- data.frame(color,size)

afc.1 <- afc(burt(x))
afc.2 <- afc(matlogic(x))

plotAll(afc.1)
plotAll(afc.2)
}
}

\keyword{multivariate}

\author{Antoine Lucas}



