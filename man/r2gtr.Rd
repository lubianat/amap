\name{r2gtr}
\alias{r2gtr}
\alias{r2atr}
\alias{r2cdt}
\title{Write to gtr, atr, cdt file format}
\description{Write data frame and hclust object to gtr atr,
  cdt files (Xcluster or Cluster output).
  Visualisation of cluster can be
  done with tools like treeview}
\usage{
r2gtr(hr,file="cluster.gtr",distance=hr$dist.method,dec='.',digits=5)
r2atr(hc,file="cluster.atr",distance=hc$dist.method,dec='.',digits=5)
r2cdt(hr,hc,data,labels=FALSE,description=FALSE,file="cluster.cdt",dec='.')
}

\details{Function \code{\link{hclust2treeview}} compute hierarchical
  clutering and export to all files at once.}

\arguments{
  \item{file}{the path of the file}
  \item{data}{a matrix (or data frame) which provides the data to put into the file} 
  \item{hr,hc}{objects of class hclust (rows and columns)}
  \item{distance}{The distance measure used. This must be one of
          `"euclidean"', `"maximum"', `"manhattan"', `"canberra"' or
          `"binary"'. Any unambiguous substring can be given.} 
  \item{digits}{number digits for precision} 
  \item{labels}{a logical value indicating whether we use the frist column as labels (NAME column for cluster file)}
  \item{description}{a logical value indicating whether we use the second column as description (DESCRIPTION column for cluster file)}
  \item{dec}{the character used in the file for decimal points}

}

\examples{
#    Create data
.Random.seed <- c(1,  416884367 ,1051235439)
m <- matrix(rep(1,3*24),ncol=3)  
m[9:16,3] <- 3 ; m[17:24,] <- 3    #create 3 groups
m <- m+rnorm(24*3,0,0.5)           #add noise
m <- floor(10*m)/10                #just one digits

# use library stats
# Cluster columns
hc <- hclust(dist(t(m)))
# Cluster rows
hr <- hclust(dist(m))

# Export files
r2atr(hc,file="cluster.atr")
r2gtr(hr,file="cluster.gtr")
r2cdt(hr,hc,m ,file="cluster.cdt")

}
\keyword{cluster}


\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}

\seealso{\code{\link{r2xcluster}}, \code{\link{xcluster2r}},\code{\link[stats]{hclust}},\code{\link{hcluster}},\code{\link{hclust2treeview}}}
