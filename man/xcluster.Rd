\name{xcluster}
\alias{xcluster}
\title{Hierarchical clustering}
\description{Performs a hierarchical cluster analysis on a set of
     dissimilarities (this function launch an external program: Xcluster).}
\usage{
xcluster(data,distance="euclidean",clean=FALSE,tmp.in="tmp.txt",tmp.out="tmp.gtr")
}

\arguments{
  \item{data}{a matrix (or data frame) which provides the data to analyze} 
  \item{distance}{The distance measure used with \emph{Xcluster}. This must be one of
    \code{"euclidean"}, \code{"pearson"} or \code{"notcenteredpearson"}.
    Any unambiguous substring can be given.}
    \item{clean}{a logical value indicating whether you want the true
    distances (\code{clean=FALSE}), or you want a clean dendrogram}
    \item{tmp.in, tmp.out}{temporary files for Xcluster}
}

\details{
  Available distance measures are (written for two vectors \eqn{x} and
  \eqn{y}): 
  \itemize{
    \item Euclidean: Usual square distance between the two vectors (2 norm).
    \item Pearson: \eqn{1- \mbox{cor}(x,y)}{1 - cor(x,y)}
    \item Pearson not centered: \eqn{1 - \frac{\sum_i x_i y_i}{\left(\sum_i x_i^2
     \sum_i y_i^2 \right) ^{1/2}}   }{1 - [ sum x_i y_i ] / sqrt[ sum x_i^2 * sum y_i^2 ] }
    }
    
Xcluster does not use usual agglomerative methods (single, average, complete),
but    compute  the distance between each 
    groups' barycenter  for the distance between two groups.

This have a problem for this kind of data:

\tabular{rll}{
A     \tab  0   \tab     0\cr
B     \tab  0   \tab     1\cr
C     \tab   0.9   \tab  0.5\cr
}

Ie: a triangular in {\bf R}$^2$, the distance between A and B is larger
than the distance between the group A,B and C (with euclidean distance).

For that case it can be useful to use \code{clean=TRUE} and that mean
that you must not consider A and B as a group without C.
}


\value{
  An object of class \bold{hclust} which describes the
  tree produced by the clustering process.
  The object is a list with components:

  \item{merge}{an \eqn{n-1} by 2 matrix.
    Row \eqn{i} of \code{merge} describes the merging of clusters
    at step \eqn{i} of the clustering.
    If an element \eqn{j} in the row is negative,
    then observation \eqn{-j} was merged at this stage.
    If \eqn{j} is positive then the merge
    was with the cluster formed at the (earlier) stage \eqn{j}
    of the algorithm.
    Thus negative entries in \code{merge} indicate agglomerations
    of singletons, and positive entries indicate agglomerations
    of non-singletons.}

  \item{height}{a set of \eqn{n-1} non-decreasing real values.
    The clustering \emph{height}: that is, the value of
    the criterion associated with the clustering
    \code{method} for the particular agglomeration.}

  \item{order}{a vector giving the permutation of the original
    observations suitable for plotting, in the sense that a cluster
    plot using this ordering and matrix \code{merge} will not have
    crossings of the branches.}

  \item{labels}{labels for each of the objects being clustered.}

  \item{call}{the call which produced the result.}

  \item{method}{the cluster method that has been used.}

  \item{dist.method}{the distance that has been used to create \code{d}
    (only returned if the distance object has a \code{"method"}
    attribute).}
}

\note{
\emph{Xcluster} is a C program made by \emph{Gavin Sherlock} that performs
hierarchical clustering, K-means and SOM. 

\emph{Xcluster} is copyrighted. 
To get or have information
about
\emph{Xcluster}: \url{http://genome-www.stanford.edu/~sherlock/cluster.html}
}

\examples{
#    Create data
.Random.seed <- c(1,  416884367 ,1051235439)
m <- matrix(rep(1,3*24),ncol=3)  
m[9:16,3] <- 3 ; m[17:24,] <- 3    #create 3 groups
m <- m+rnorm(24*3,0,0.5)           #add noise
m <- floor(10*m)/10                #just one digits


# And once you have Xcluster program:
#
#h <- xcluster(m)
#
#plot(h) 
}
\keyword{cluster}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}
\seealso{\code{\link{r2xcluster}}, \code{\link{xcluster2r}},\code{\link[stats]{hclust}}, \code{\link{hcluster}}}
