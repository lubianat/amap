\name{plot.acp}
\alias{plot.acp}
\alias{biplot.acp}
\alias{plot2.acp}

\title{Graphics for Principal component Analysis}
\description{Graphics for Principal component Analysis / Graphiques pour l'analyse en composantes principales, "comme Carlier"}
\usage{
plot.acp(x,i=1,j=2,text=TRUE,label='Composante ',col='darkblue',
main='ACP des individus',...)
biplot.acp(x,i=1,j=2,label='Composante ',col='darkblue',length=0.1,
main='ACP des variables',...)
plot2.acp(x,pourcent=FALSE,eigen=TRUE,label='Comp.',col='lightgrey',
main='Eboulis des valeurs propres',ylab='Valeurs propres')}

\arguments{
\item{x}{Result of acp / princomp}
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
\item{\dots}{cex, pch, and other options; see points.}
}




\value{
   Graphics: 

\code{plot.acp} PCA for lines / ACP des individus

\code{plot.acp} PCA for columns / ACP des variables

\code{plot2.acp} Eigen values diagram / Eboulis des valeurs propres
}



\examples{
data(lubisch)
lubisch <- lubisch[,-c(1,8)]
p <- acp(lubisch)
plot(p)
}
\keyword{multivariate}

\author{Antoine Lucas, \url{http://genopole.toulouse.inra.fr/~lucas/amap}}
\seealso{\link{acp},\link{acpgen}, \link[mva]{princomp}}
