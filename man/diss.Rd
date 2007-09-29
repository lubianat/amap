\name{diss}
\alias{diss}

\title{Compute a dissimilarity matrix}
\description{
  Compute a dissimilarity matrix from a data set (containing only factors).
  
}
\usage{
diss(x, w=rep(1,ncol(x)) )
}

\arguments{
  \item{x}{A matrix or data frame containing only factors.}
  \item{w}{A vector of weight, by default each variable has got same weight}
}
\details{
%%  Dissimilarity between individual i and j are number of common factors
%%  substracted to number of different factors.
   Case of N individuals described by P categorical variables:
   each element (i,j) of the signed similarities array is computed
   by sommation over the P variables of the contributions of each
   variable, multiplied by the weight of the variable. The contribution
   of a given categorical variable is +1 if the individual i and j
   are in the same class, and is -1 if they are not.
}
\value{
  A dissimilarity matrix.
}

\author{Antoine Lucas}

\seealso{ \code{\link{Dist}}, \code{\link{pop}} }
\examples{
data <-
matrix(c(1,1,1,1,1
        ,1,2,1,2,1
        ,2,3,2,3,2
        ,2,4,3,3,2
        ,1,2,4,2,1
        ,2,3,2,3,1),ncol=5,byrow=TRUE)

diss(data)

## With weights
diss(data,w=c(1,1,2,2,3))

}
\keyword{multivariate}
\keyword{cluster}

