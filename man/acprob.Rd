\name{acprob}
\alias{acprob}
\title{Robust principal component analysis}
\description{Robust principal component analysis / Analyse en composantes principales robuste}

\usage{
acprob(x,h,center=TRUE,reduce=TRUE,kernel="gaussien")
}

\arguments{
\item{x}{Matrix  / data frame}
\item{h}{Scalar: bandwidth of the Kernel}
\item{kernel}{The kernel used. This must be one of '"gaussien"', 
'"quartic"', '"triweight"', '"epanechikov"' , 
'"cosinus"' or '"uniform"' } 
\item{center}{A logical value indicating whether we center data}
\item{reduce}{A logical value indicating whether we "reduce" data i.e. 
divide each column by standard deviation}
}


\details{
  \code{acpgen}  compute robust pca. i.e. spectral analysis of a robust
  variance instead of usual variance. Robust variance: see
  \code{\link[amap]{varrob}}
  }

  \value{
  An object of class \bold{acp} 
  The object is a list with components:

  \item{sdev}{the standard deviations of the principal components.}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix
    whose columns contain the eigenvectors).  This is of class
    \code{"loadings"}: see \code{\link[mva]{loadings}} for its \code{print}
    method.}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied
    data on the principal components.}
  \item{eig}{Eigen values}
}
\keyword{multivariate}

\author{Antoine Lucas, \url{http://genopole.toulouse.inra.fr/~lucas/amap}}
\seealso{\link{acp} \link[mva]{princomp} }
