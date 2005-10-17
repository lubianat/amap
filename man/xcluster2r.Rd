\name{xcluster2r}
\alias{xcluster2r}
\title{Importing Xcluster/Cluster output}
\description{Converting Xcluster/Cluster output (.gtr or .atr) to R hclust file}
\usage{
xcluster2r(file,distance="euclidean",labels=FALSE,fast=FALSE,clean=FALSE,
           dec='.')
}

\arguments{
  \item{file}{the path of a \emph{Xcluster/Cluster} file (.gtr or .atr)}
  \item{distance}{The distance measure used with \emph{Xcluster/Cluster}. This must be one of
    \code{"euclidean"}, \code{"pearson"} or \code{"notcenteredpearson"}.
    Any unambiguous substring can be given.}
  \item{labels}{a logical value indicating whether we use 
    labels values (in the .cdt file) or not.}
  \item{fast}{a logical value indicating whether we reorganize data
    like R (\code{Fast=FALSE}) or we let them like \emph{Xcluster/Cluster} did} 
  \item{clean}{a logical value indicating whether you want the true
    distances (\code{clean=FALSE}), or you want a clean dendrogram (see
  details below).}
  \item{dec}{the character used in the file for decimal points}
}

\details{
See \bold{xcluster} for more details.
}


\value{
  An object of class \bold{hclust} which describes the
  tree produced by the clustering process.
}

\note{
\emph{Xcluster} is a C program made by \emph{Gavin Sherlock} that performs
hierarchical clustering, K-means and SOM. 

\emph{Xcluster} is copyrighted. 
To get or have information
about
\emph{Xcluster}: \url{http://genome-www.stanford.edu/~sherlock/cluster.html}

\emph{Cluster} is a program made by \emph{Michael Eisen} that performs
hierarchical clustering, K-means and SOM. 

\emph{Cluster} is copyrighted. 
To get or have information
about
\emph{Cluster}: \url{http://rana.lbl.gov/EisenSoftware.htm}
}

\examples{
#    Create data
.Random.seed <- c(1,  416884367 ,1051235439)
m <- matrix(rep(1,3*24),ncol=3)  
m[9:16,3] <- 3 ; m[17:24,] <- 3    #create 3 groups
m <- m+rnorm(24*3,0,0.5)           #add noise
m <- floor(10*m)/10                #just one digits

r2xcluster(m)


# And once you have Xcluster program:
#
#system('Xcluster -f xcluster.txt -e 0 -p 0 -s 0 -l 0')
#h <- xcluster2r('xcluster.gtr')
#plot(h,hang=-1)
}
\keyword{file}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}
\seealso{\code{\link{xcluster}}, \code{\link{r2xcluster}}, \code{\link[stats]{hclust}}, \code{\link{hcluster}}}
