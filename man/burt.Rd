\name{burt}
\alias{burt}
\alias{matlogic}
\title{Compute burt table from a factor dataframe.}

\usage{
burt(x)
matlogic(x)
}

\arguments{
  \item{x}{A dataframe that contents only factors}
}

\description{
  %matlogic: tableau disjonctif complet
  %burt
  matlogic returns for all variables a matrix of logical values
  for each levels. burt is defined as t(matlogic).matlogic
}

\examples{
color <- as.factor(c('blue','red','red','blue','red'))
size <- as.factor(c('large','large','small','medium','large'))
x <- data.frame(color,size)

matlogic(x)
##  color.blue color.red size.large size.medium size.small
##1          1         0          1           0          0
##2          0         1          1           0          0
##3          0         1          0           0          1
##4          1         0          0           1          0
##5          0         1          1           0          0

burt(x)
##              color.blue color.red size.large size.medium size.small
##  color.blue           2         0          1           1          0
##  color.red            0         3          2           0          1
##  size.large           1         2          3           0          0
##  size.medium          1         0          0           1          0
##  size.small           0         1          0           0          1
}

\keyword{multivariate}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}



