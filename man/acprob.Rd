\name{acprob}
\alias{acprob}
\title{Robust principal component analysis}
\description{Robust principal component analysis}

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
  \code{\link{varrob}}
  }

  \value{
  An object of class \bold{acp} 
  The object is a list with components:

  \item{sdev}{the standard deviations of the principal components.}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix
    whose columns contain the eigenvectors).  This is of class
    \code{"loadings"}: see \code{\link[stats]{loadings}} for its \code{print}
    method.}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied
    data on the principal components.}
  \item{eig}{Eigen values}
}
\keyword{multivariate}

\references{
  H. Caussinus, M. Fekri, S. Hakam and A. Ruiz-Gazen, \emph{A monitoring
  display of multivariate outliers} Computational Statistics & Data
  Analysis,
  Volume 44, Issues 1-2, 28 October 2003, Pages 237-252

  Caussinus, H and Ruiz-Gazen, A. (1993): \emph{Projection Pursuit and
  Generalized Principal Component Analyses, in New Directions in
  Statistical Data Analysis and Robustness} (eds.  Morgenthaler et
  al.), pp. 35-46. Birkäuser Verlag Basel.

  Caussinus, H. and Ruiz-Gazen, A. (1995). \emph{Metrics for Finding Typical
  Structures by Means of Principal Component Analysis. In Data Science
  and its Applications} (eds Y. Escoufier and C. Hayashi),
  pp. 177-192. Tokyo: Academic Press.
}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}
\seealso{\link[stats]{princomp} \link{acpgen}}
