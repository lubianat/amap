#-------------------------------------------------------
#
#  Created       : 29/10/02
#  Last Modified : Time-stamp: <2002-11-12 09:59:15 lucas>
#
#  Description   : Principal component analysis
#                  
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#
#  Licence       : GPL 
#
#-------------------------------------------------------



acp <- function(x,center=TRUE,reduce=TRUE)
{   
    x    <- as.matrix(x)
    x    <- scale(x ,center = center, scale = FALSE)
    if (reduce == TRUE)
         {
          x    <- apply(x,2,function(u) { u/sd(u)}) 
         }
    EIG  <- eigen( t(x) %*% x,symmetric=TRUE) 
    V    <- EIG$vector    # ou bien: V=svd(x)$v
    val  <- sqrt(EIG$values)

    scores <- x %*% V
    sdev   <- apply(scores,2,sd)    
    res  <- list(eig=val,sdev=sdev,scores=scores,loadings=V)
    class(res) <- "acp"
    res
}
pca <- acp

print.acp <- function(x, ...)
{
    #cat("Call:\n"); dput(x$call)
    cat("\nStandard deviations:\n")
    print(x$sdev, ...)
    cat("\nEigen values:\n")
    print(x$eig, ...)
    invisible(x)
}


# 
#   SECTION GRAPHIQUES
#

plot.acp <- function(x,i=1,j=2,text=TRUE,label='Composante ',col='darkblue',main='ACP des individus',...)
{
    U    <- x$scores
    XLAB <- paste(label,i)
    YLAB <- paste(label,j)
    plot(U[,i],U[,j],col='white',xlab=XLAB,ylab=YLAB,main=main)
    if(text){
        text(U[,i],U[,j],col=col,...)   
    }
    else{
        points(U[,i],U[,j],col=col,...) 
    }
}

biplot.acp <- function(x,i=1,j=2,label='Composante ',col='darkblue',length=0.1,main='ACP des variables')
{
    U    <- x$loadings
    LIM  <- c(-1.3,1.3)
    XLAB <- paste(label,i)
    YLAB <- paste(label,j)
    # PLOT DES AXES
    plot.default(0,col='white',xlim=LIM,ylim=LIM,xlab=XLAB,ylab=YLAB,main=main)

    # PLOT DU NOM DES FLECHES
    text(U[,i]*1.3,U[,j]*1.3,labels=dimnames(U)[[1]],col=col)   

    # PLOT DES FLECHES
    arrows(0,0,U[,i],U[,j],length = length,col=col)

    # CERCLE
    t2p <- 2 * pi * seq(0,1, length = 200)
    xc <- cos(t2p)
    yc <- sin(t2p)
    lines(xc,yc,col='darkblue')
}

# Graphique: Eboulis des valeurs propres
plot2.acp <- function(x,pourcent=FALSE,eigen=TRUE,label='Comp.',col='lightgrey',main='Eboulis des valeurs propres',ylab='Valeurs propres')
{
    if(eigen){ U <- x$eig }
    else { U <- x$sdev }

    if(pourcent){U <- U/sum(U) }
    n     <- length(U)
    names <- paste(label,1:n)
    barplot(U,main=main,ylab=ylab,col=col,names.arg=names)
}


