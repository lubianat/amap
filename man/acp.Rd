\name{acp}
\alias{acp}
\alias{pca}
\alias{print.acp}

\title{Principal component analysis}
\description{Principal component analysis / Analyse en composantes principales}
\usage{
acp(x,center=TRUE,reduce=TRUE)
pca(x,center=TRUE,reduce=TRUE)
print.acp(x, ...) 
}

\arguments{
\item{x}{Matrix  / data frame}
\item{center}{a logical value indicating whether we center data}
\item{reduce}{a logical value indicating whether we "reduce" data i.e. 
divide each column by standard deviation}
\item{\dots}{arguments to be passed to or from other methods.}
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

\references{
  A. Carlier
  Analyse des données Multidimensionnelles
  \url{http://www.lsp.ups-tlse.fr/Carlier/enseignement.html}
}


\examples{
data(lubisch)
lubisch <- lubisch[,-c(1,8)]
p <- acp(lubisch)
plot(p)
}
\keyword{multivariate}

\author{Antoine Lucas, \url{http://genopole.toulouse.inra.fr/~lucas/amap}}
\seealso{\link{plot.acp},\link{acpgen}, \link[mva]{princomp} }
