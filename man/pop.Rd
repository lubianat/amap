\name{pop}
\alias{pop}
\title{Optimal Partition (classification).}

\usage{
pop(x,fmbvr=TRUE,triabs=TRUE,allsol=TRUE)
}

\arguments{
  \item{x}{A dissimilarity matrix}
  \item{fmbvr}{Logical, TRUE: look for the exact solution}
  \item{triabs}{Logical, TRUE: try to init with absolute values}
  \item{allsol}{Logical, TRUE all solutions, FALSE only one solution}
}

\description{
 Classification: Computing an Optimal Partition
from Weighted Categorical Variables
or from an Array of Signed Similarities.  
}

\references{
 Theory is explained at \url{http://petitjeanmichel.free.fr/itoweb.petitjean.class.html}
 
 Marcotorchino F.
 \emph{Agr\'egation des similarit\'es en classification automatique.}
 Th\'ese de Doctorat d'Etat en Math\'ematiques,
 Universit\'e Paris VI, 25 June 1981.

 Petitjean M.
 \emph{Agr\'egation des similarit\'es: une solution oubli\'ee.}
 RAIRO Oper. Res. 2002,36[1],101-108.

  }

  

\examples{

## pop from a data matrix
data <-
matrix(c(1,1,1,1,1
        ,1,2,1,2,1
        ,2,3,2,3,2
        ,2,4,3,3,2
        ,1,2,4,2,1
        ,2,3,2,3,1),ncol=5,byrow=TRUE)



pop(diss(data))


## pop from a dissimilarity matrix

d <-2 * matrix(c(9,  8,  5,  7,  7,  2
,  8,  9,  2,  5,  1,  7
,  5,  2,  9,  8,  7,  1
,  7,  5,  8,  9,  3,  2
,  7,  1,  7,  3,  9,  6
,  2,  7,  1,  2,  6,  9),ncol=6,byrow=TRUE) - 9

pop(d)


\dontrun{
d <- 2 * matrix(c(57, 15, 11, 32,  1, 34,  4,  6, 17,  7
, 15, 57, 27, 35, 27, 27, 20, 24, 30, 15
, 11, 27, 57, 25, 25, 20, 34, 25, 17, 15
, 32, 35, 25, 57, 22, 44, 13, 22, 30, 11
,  1, 27, 25, 22, 57, 21, 28, 43, 20, 13
, 34, 27, 20, 44, 21, 57, 18, 27, 21,  8
,  4, 20, 34, 13, 28, 18, 57, 31, 28, 13
,  6, 24, 25, 22, 43, 27, 31, 57, 30, 15
, 17, 30, 17, 30, 20, 21, 28, 30, 57, 12
,  7, 15, 15, 11, 13,  8, 13, 15, 12, 57),ncol=10,byrow=TRUE) - 57

pop(d)


}

}

\keyword{multivariate}


\author{Michel Petitjean,
  \url{http://petitjeanmichel.free.fr/itoweb.petitjean.class.html}
  
  R port by Antoine Lucas,
}



