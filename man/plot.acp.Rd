\name{plot.acp}
\alias{plot.acp}
\alias{biplot.acp}
\alias{plot2}
\alias{plotAll}

\title{Graphics for Principal component Analysis}
\description{Graphics for Principal component Analysis}
\usage{
plot.acp(x,i=1,j=2,text=TRUE,label='Composants',col='darkblue',
main='Individuals PCA',variables=TRUE,...)
biplot.acp(x,i=1,j=2,label='Composants',col='darkblue',length=0.1,
main='Variables PCA',circle=TRUE,...)
plot2(x,pourcent=FALSE,eigen=TRUE,label='Comp.',col='lightgrey',
main='Scree Graph',ylab='Eigen Values')
plotAll(x)
}

\arguments{
\item{x}{Result of acp or princomp}
\item{i}{X axis}
\item{j}{Y axis}
\item{text}{a logical value indicating whether we use text or points for plot}
\item{pourcent}{a logical value indicating whether we use pourcentage of 
values}
\item{eigen}{a logical value indicating whether we use eigen values or
standard deviation}
\item{label}{label for X and Y axis}
\item{col}{Color of plot}
\item{main}{Title of graphic}
\item{ylab}{Y label}
\item{length}{length of arrows}
\item{variables,circle}{a logical value indicating whether we display
  circle or variables}
\item{\dots}{cex, pch, and other options; see points.}
}




\value{
   Graphics: 

\code{plot.acp} PCA for lines (individuals)

\code{plot.acp} PCA for columns (variables)

\code{plot2} Eigen values diagram (Scree Graph)

\code{plotAll} Plot both 3 graphs
}



\examples{
data(lubisch)
lubisch <- lubisch[,-c(1,8)]
p <- acp(lubisch)
plotAll(p)
}
\keyword{multivariate}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}
\seealso{\link{acpgen},\link{acprob}, \link[stats]{princomp}}
