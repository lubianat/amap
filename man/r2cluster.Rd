\name{r2cluster}
\alias{r2cluster}
\title{Write to Cluster file format}
\description{Converting data to Cluster format}
\usage{
r2cluster(data,labels=FALSE,colname="ACC",description=FALSE,
          file="cluster.txt",dec='.')
}

\arguments{
  \item{file}{the path of the file}
  \item{data}{a matrix (or data frame) which provides the data to put into the file} 
  \item{labels}{a logical value indicating whether we use the frist
    column as labels (ACC column in cluster file)}
  \item{colname}{a character string
    indicating  what kind of objects are in each row.
    YORF, MCLID, CLID, ACC can be used: see details.
  }
   \item{description}{a logical value indicating whether we use the second column as description (NAME column for cluster file)}
   \item{dec}{the character used in the file for decimal points}  
}

\details{
  Software \emph{Cluster}, made by \emph{M. Eisen} needs formatted
  input data like:

  \preformatted{
ACC     NAME    GWEIGHT GORDER  V3      V4      V5
EWEIGHT                         1       1       1
gbk01   Gene1   1       1       0.9     0.4     1.4
gbk02   Gene2   1       2       0.6     0.2     0.2
gbk03   Gene3   1       3       1.6     1.1     0.9
gbk04   Gene4   1       4       0.4     1       1
  }

  First field of first line (i.e "ACC") is a special field, that tells
  program what kind of objects are in each row.

  
  Four special values are defined with web link (when visualize with
  TreeView):
  \itemize{
    \item  YORF
    \url{http://genome-www.stanford.edu/cgi-bin/dbrun/SacchDB?find+Locus+\%22UNIQID\%22}
    \item MCLID
    \url{http://genome.rtc.riken.go.jp/cgi-bin/getseq?g+R+UNIQID}
    
    \item CLID
    \url{http://genome-www4.stanford.edu/cgi-bin/SMD/source/sourceRes\-ult?op\-tion=CloneID&criteria1=IMAGE:UNIQID&choice=cDNA}
    
    \item ACC
    \url{http://genome-www4.stanford.edu/cgi-bin/SMD/source/sourceRes\-ult?op\-tion=Number&criteria=UNIQID&choice=cDNA}
  }

  Line begining with \code{EWEIGHT} gives weights for each
  column (variable). Column \code{GWEIGHT} gives weights for each line
  (individuals).

}

\note{
\emph{Cluster} is a program made by \emph{M. Eisen} that performs
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

r2cluster(m)

}
\keyword{file}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}

\seealso{\code{\link{xcluster}}, \code{\link{r2xcluster}}, \code{\link[stats]{hclust}}}
