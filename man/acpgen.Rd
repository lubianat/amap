\name{acpgen}
\alias{acpgen}
\alias{K}
\alias{W}
\alias{WsansC}


\title{Generalised principal component analysis}
\description{Generalised principal component analysis / Analyse en composantes principales généralisée}

\usage{
acpgen(x,h1,h2,center=TRUE,reduce=TRUE,kernel="gaussien")
K(u,kernel="gaussien")
W(x,h,D=NULL,kernel="gaussien")
WsansC(x,h,D=NULL,kernel="gaussien")
}

\arguments{
\item{x}{Matrix  / data frame}
\item{h}{Scalar: bandwidth of the Kernel}
\item{h1}{Scalar: bandwidth of the Kernel for W}
\item{h2}{Scalar: bandwidth of the Kernel for U}
\item{kernel}{The kernel used. This must be one of '"gaussien"', 
'"quartic"', '"triweight"', '"epanechikov"' , 
'"cosinus"' or '"uniform"' } 
\item{center}{A logical value indicating whether we center data}
\item{reduce}{A logical value indicating whether we "reduce" data i.e. 
divide each column by standard deviation}
\item{D}{A product scalar matrix / une matrice de produit scalaire}
\item{u}{Vector}
}

\details{
\code{acpgen} compute generalised pca. i.e. spectral analysis of 
\eqn{U_n . W_n^{-1}}{Un / Wn}, and project \eqn{X_i}{Xi} with
\eqn{W_n^{-1}}{1/Wn} on the principal vector sub-spaces.

\eqn{X_i}{Xi} a column vector of \eqn{p} variables of individu  \eqn{i} 
(input data)

\code{W} compute estimation of noise in the variance. 
\deqn{W_n=\frac{\sum_{i=1}^{n-1}\sum_{j=i+1}^{n}K(||X_i-X_j||_{V_n^{-1}}/h)(X_i-X_j)(X_i-X_j)'}{\sum_{i=1}^{n-1}\sum_{j=i+1}^{n}K(||X_i-X_j||_{V_n^{-1}}/h)}}{W: see latex doc}


with \eqn{V_n}{Vn} variance estimation;

\code{U} compute robust variance. \eqn{U_n^{-1} = S_n^{-1} - 1/h V_n^{-1}}{1/Un
 = 1/Sn - 1 / (h Vn)}


\deqn{S_n=\frac{\sum_{i=1}^{n}K(||X_i||_{V_n^{-1}}/h)(X_i-\mu_n)(X_i-\mu_n)'}{\sum_{i=1}^nK(||X_i||_{V_n^{-1}}/h)}}{S: see latex doc}

with \eqn{\mu_n} estimator of the mean.

%\eqn{\delta=K(|| X_i-X_j||)} with K: a kernel.

\code{K} compute   kernel, i.e.

gaussien: \deqn{\frac{1}{\sqrt{2\pi}}  e^{-u^2/2}}{
               1/sqrt(2pi)  exp(-u^2/2)}

quartic: \deqn{\frac{15}{16}(1-u^2)^2 I_{|u|\leq 1}  }{
              15/16 (1-u^2)^2  if |u| < 1}

triweight: \deqn{\frac{35}{32}(1-u^2)^3 I_{|u|\leq 1}  }{
              35/32 (1-u^2)^3  if |u| < 1}

epanechikov: \deqn{\frac{3}{4}(1-u^2) I_{|u|\leq 1}  }{
              3/4 (1-u^2)  if |u| < 1}


cosinus: \deqn{\frac{\pi}{4}\cos(\frac{\pi}{2}u) I_{|u|\leq 1}  }{
               pi/4 cos (u * pi/2)  if |u| < 1}

}

\value{
  An object of class \bold{acp} 
  The object is a list with components:

  \item{sdev}{the standard deviations of the principal components.}
  \item{loadings}{the matrix of variable loadings (i.e., a matrix
    whose columns contain the eigenvectors).  This is of class
    \code{"loadings"}: see \code{\link[stats]{loadings}} for its \code{print}
    method.}
  \item{scores}{if \code{scores = TRUE}, the scores of the supplied
    data on the principal components.}
  \item{eig}{Eigen values}
}


\examples{
data(lubisch)
lubisch <- lubisch[,-c(1,8)]
p <- acpgen(lubisch,h1=1,h2=1/sqrt(2))
plot(p,main='ACP robuste des individus')

# See difference with acp

p <- acp(lubisch)
plot(p)

# Other Data
.Random.seed <- c(1,  416884367 ,1051235439)
n <- 100   # Number of observations
m <- 0     # Mean of 1 group
temp <- cbind(rnorm(n,m,1),rnorm(n,0,10))
m <- 10
mat <- rbind(temp,cbind(rnorm(n,m,1),rnorm(n,0,10)))
mat <- scale(mat,TRUE,TRUE)
rm (m,n,temp) 

}

\references{
  H. Caussinus, M. Fekri, S. Hakam and A. Ruiz-Gazen, \emph{A monitoring
  display of multivariate outliers} Computational Statistics & Data
  Analysis,
  Volume 44, Issues 1-2, 28 October 2003, Pages 237-252

  Caussinus, H and Ruiz-Gazen, A. (1993): \emph{Projection Pursuit and
  Generalized Principal Component Analyses, in New Directions in
  Statistical Data Analysis and Robustness} (eds.  Morgenthaler et
  al.), pp. 35-46. Birkäuser Verlag Basel.

  Caussinus, H. and Ruiz-Gazen, A. (1995). \emph{Metrics for Finding Typical
  Structures by Means of Principal Component Analysis. In Data Science
  and its Applications} (eds Y. Escoufier and C. Hayashi),
  pp. 177-192. Tokyo: Academic Press.
}


\keyword{multivariate}

\author{Antoine Lucas, \url{http://mulcyber.toulouse.inra.fr/projects/amap/}}
\seealso{\link{acp} \link{acprob} \link[stats]{princomp} }



