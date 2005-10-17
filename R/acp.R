#-------------------------------------------------------
#
#  Created       : 29/10/02
#  Last Modified : Time-stamp: <2005-10-02 17:55:16 antoine>
#
#  Description   : Principal component analysis
#                  
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#
#  Licence       : GPL 
#
#-------------------------------------------------------




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

plot.acp <- function(x,i=1,j=2,text=TRUE,label='Composants',col='darkblue',main='Individuals PCA',...)
{
    U    <- x$scores
    XLAB <- paste(label,i)
    YLAB <- paste(label,j)
    plot.new()
    plot.window(range(U[,i]),range(U[,j]))
    axis(1,label=TRUE,tick=TRUE)
    axis(2,label=TRUE,tick=TRUE)
    box()
    title(xlab=XLAB,ylab=YLAB,main=main)
    if(text){
        text(U[,i],U[,j],col=col,...)   
    }
    else{
        points(U[,i],U[,j],col=col,...) 
    }
}

biplot.acp <- function(x,i=1,j=2,label='Composants',col='darkblue',length=0.1,main='Variables PCA',...)
{
    U    <- x$loadings
    LIM  <- c(-1.3,1.3)
    XLAB <- paste(label,i)
    YLAB <- paste(label,j)

    # PLOT DES AXES
    plot.new()
    plot.window(LIM,LIM)
    axis(1,label=TRUE,tick=TRUE)
    axis(2,label=TRUE,tick=TRUE)
    box()
    title(xlab=XLAB,ylab=YLAB,main=main)


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
plot2.acp <- function(x,pourcent=FALSE,eigen=TRUE,label='Comp.',col='lightgrey',main='Scree Graph',ylab='Eigen Values')
{
    if(eigen){ U <- x$eig }
    else { U <- x$sdev }

    if(pourcent){U <- U/sum(U) }
    n     <- length(U)
    names <- paste(label,1:n)
    barplot(U,main=main,ylab=ylab,col=col,names.arg=names)
}


