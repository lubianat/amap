#-------------------------------------------------------
#
#  Created       : 30/10/02
#  Last Modified : Time-stamp: <2002-11-25 10:19:04 lucas>
#
#  Description   : Robust principal component analysis
#                  
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#
#  Licence       : GPL 
#
#-------------------------------------------------------

K <- function(u,kernel="gaussien") {
    switch(kernel,
        gaussien = (2*pi)^(-1/2) * exp(-u^2/2),
        quartic   = 15/16 * (1-u^2)^2 * (abs(u)<1),
        triweight = 35/32 * (1-u^2)^3 * (abs(u)<1),
        epanechikov = 3/4 * (1-u^2) *   (abs(u)<1),
        cosinus = pi/4 * cos (u*pi/2) * (abs(u)<1),
        uniform = 1/2 * (abs(u)<1),
    )
}

# Variance locale
W <- function(x,h,D=NULL,kernel="gaussien")
{
    x   <- as.matrix(x)
    n   <- dim(x)[1]
    p   <- dim(x)[2]
    if (is.null(D)) {
        D <- diag(1,p)
    }
    x <- as.vector(x)
    D <- as.vector(D)
    kernel <- substr(kernel,1,1)

    matrix(.C("W",as.double(x),as.double(h),as.double(D),as.integer(n),as.integer(p),as.character(kernel),res=double(p*p))$res,p)
}


WsansC <- function(x,h,D=NULL,kernel="gaussien")
{
    x   <- as.matrix(x)
    n   <- dim(x)[1]
    p   <- dim(x)[2]
    if (is.null(D)) {
        D <- diag(1,p)
    }

    som    <- 0
    result <- 0
    for(i in 1:(n-1) )
        {
        for(j in (i+1):n)
            {
            Delta <- x[i,]-x[j,]
            norm <- sqrt(t(Delta) %*%D %*%Delta)
            k <- K ( norm /h,kernel)   # K ( |Delta|/h )
            k <- as.numeric(k)
            som <- som + k
            result <- result + k * Delta %*% t(Delta)
            }
        }
    result /   som
}

U <- function(x,h,D=NULL,kernel="gaussien")
{
    x   <- as.matrix(x)
    x   <- scale(x, center = TRUE, scale = FALSE)
    n   <- dim(x)[1]
    p   <- dim(x)[2]
    if (is.null(D)) {
        D <- diag(1,p)
    }
    x <- as.vector(x)
    D <- as.vector(D)
    kernel <- substr(kernel,1,1)

    S <- matrix(.C("U",as.double(x),as.double(h),as.double(D),as.integer(n),as.integer(p),as.character(kernel),res=double(p*p))$res,p)
    Sinv <- solve(S)
    solve ( Sinv - D / h)
}


UsansC <- function(x,h,D=NULL,kernel="gaussien")
{
    n   <- dim(x)[1]
    p   <- dim(x)[2]
    if (is.null(D)) {
        D <- diag(1,p)
    }
    x   <- as.matrix(x)
    x   <- scale(x ,center = TRUE, scale = FALSE)
    som <- 0
    res <- 0
    for (i in 1:n )
    {
        k <- K (  sqrt(t(x[i,]) %*%D %*% x[i,]) /h,kernel)
        k <- as.numeric(k)
        res <- res + k * x[i,] %*% t( x[i,])
        som <- som + k
    }
    S <- res / som
    Sinv <- solve(S)
    solve ( Sinv -  D / h )

}


acpgen <- function(x,h1,h2,center=TRUE,reduce=TRUE,kernel="gaussien")
{
    # CENTRONS, ET REDUISONS
    x    <- as.matrix(x)
    x    <- scale(x ,center = center, scale = FALSE)
    if (reduce == TRUE)
         {
          x    <- apply(x,2,function(u) { u/sd(u)}) 
         }

    # ESTIMATION DE W et U
    n <- dim(x)[1]
    VarInv   <- solve(var(x)*(n-1)/n) # solve= inverser
    leU    <- U(x,h1,D=VarInv,kernel=kernel)
    leW    <- W(x,h2,D=VarInv,kernel=kernel)
    Winv   <- solve(leW) 


    # anal. spec de Var.W^-1 :
    EIG    <- eigen(leU %*% Winv)  
    V      <- EIG$vector

    #EIG    <- eigen( x %*% Winv %*% t(x)  )
    #U      <- EIG$vector
    #n      <- dim(x)[1]
    #p      <- dim(x)[2]
    #S      <- diag(Re(EIG$values),n)   
    #S1     <- diag(Re(1/EIG$values),n)
    #S      <- sqrt(S[,1:p])
    #S1     <- sqrt(S1[,1:p])
    #V      <- t(x)%*% U%*% S1
    # X=U.S.V' -> V = X' U S^-1
    

    # AFFICHAGE DES RESULTATS
    scores <- x %*% Winv %*% V
    eig    <- sqrt(EIG$values)
    sdev   <- apply(scores,2,sd)    
    res    <- list(eig=eig,sdev=sdev,scores=scores,loadings=V)
    class(res) <- "acp"
    res
}
