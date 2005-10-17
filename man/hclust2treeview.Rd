\name{hclust2treeview}
\alias{hclust2treeview}
\title{Hierarchical clustering and treeview export}
\description{
  This function compute hierachical clustering with function
  hcluster and export cluster to treeview files format.
}

\usage{
hclust2treeview(x,file="cluster.cdt",method = "euclidean",
    link = "complete",keep.hclust=FALSE)
}

\details{This function producte all 3 files needed by treeview, with
  extentions: cdt, gtr, atr.}

\arguments{
 \item{x}{numeric matrix or a data frame or an object of class
    "exprSet".}
  \item{file}{File name of export file.}
  \item{method}{the distance measure to be used. This must be one of
    \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
    \code{"canberra"} \code{"binary"} \code{"pearson"},
    \code{"correlation"} or  \code{"spearman"}.
    Any unambiguous substring can be given.}
  \item{link}{the agglomeration method to be used. This should
    be (an unambiguous abbreviation of) one of
    \code{"ward"}, \code{"single"}, \code{"complete"},
    \code{"average"}, \code{"mcquitty"}, \code{"median"} or
    \code{"centroid"}.}
  \item{keep.hclust}{if TRUE: function returns a list of 2 objects of
    class hclust}
}

\value{
  if keep.hclust = FALSE, function return 1.
  else function returns 2 objects of class hclust, first: hierachical
  clustering by rows, second: hierarchical clustering by columns
}



\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}

\seealso{\code{\link{hcluster}}}

\examples{

data(USArrests)
hclust2treeview(USArrests,file="cluster.cdt")

}


\keyword{cluster}




