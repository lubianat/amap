\name{r2xcluster}
\alias{r2xcluster}
\title{Write to Xcluster file format}
\description{Converting data to Xcluster format}
\usage{
r2xcluster(data,labels=FALSE,description=FALSE,file="xcluster.txt")
}

\arguments{
  \item{file}{the path of the file}
  \item{data}{a matrix (or data frame) which provides the data to put into the file} 
   \item{labels}{a logical value indicating whether we use the frist column as labels (NAME column for xcluster file)}
   \item{description}{a logical value indicating whether we use the second column as description (DESCRIPTION column for cluster file)}
}


\note{
\emph{Xcluster} is a C program made by \emph{Gavin Sherlock} that performs
hierarchical clustering, K-means and SOM. 

\emph{Xcluster} is copyrighted. 
To get or have information
about
\emph{Xcluster}: \url{http://genome-www.stanford.edu/~sherlock/cluster.html}
}

\details{
  Software \emph{Xcluster}, made by \emph{G. Sherlock} needs formatted
  input data like:

  \preformatted{
NAME    DESCRIPTION     GWEIGHT V2      V3      V4
EWEIGHT                         1       1       1
gbk01   Gene1           1       0.9     0.4     1.4
gbk02   Gene2           1       0.6     0.2     0.2
gbk03   Gene3           1       1.6     1.1     0.9
gbk04   Gene4           1       0.4     1       1
  }

  Line begining with \code{EWEIGHT} gives weights for each
  column (variable). Column \code{GWEIGHT} gives weights for each line
  (individuals).
  

}

\examples{
##   Create data
.Random.seed <- c(1,  416884367 ,1051235439)
m <- matrix(rep(1,3*24),ncol=3)  
m[9:16,3] <- 3 ; m[17:24,] <- 3    #create 3 groups
m <- m+rnorm(24*3,0,0.5)           #add noise
m <- floor(10*m)/10                #just one digits

r2xcluster(m)

## And once you have Xcluster program:

\dontrun{
system('Xcluster -f xcluster.txt -e 0 -p 0 -s 0 -l 0')
h <- xcluster2r('xcluster.gtr')
plot(h,hang=-1)
}
}
\keyword{file}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}

\seealso{\code{\link{xcluster}}, \code{\link{xcluster2r}}, \code{\link[stats]{hclust}}, \code{\link{hcluster}}}
