
##
matlogic <- function(x)
{
  n=nrow(x)
  m=ncol(x)
  nblev <- apply(x,2,function(u){nlevels(as.factor(u))})

  ## Keep names....
  rownames <- rownames(x)
  colnames <- colnames(x)
  i <- 0
  colnamesnew <- c(apply(x,2,function(u){ i<<- i+1;paste(colnames[i],levels(as.factor(u)),sep=".")}),recursive=TRUE)
  

  k <- sum(nblev)
  res <- as.integer(matrix(0,ncol=k,nrow=n))
  x <- c(x,recursive=TRUE)
  
  res <- .C("matind",
            as.integer(nblev),
            as.integer(x),
            res=res,
            n,
            m,
            as.integer(k),
            PACKAGE="amap")

  res <- matrix(res$res,ncol=k)
  rownames(res) <- rownames
  colnames(res) <- colnamesnew
  res
  
}


burt <- function(x)
  {
    ind <- matlogic(x)
    t(ind) %*% ind

  }



## x: table de burt, ou table
afc <- function (x)
  {
    f  <- as.matrix(x/sum(x))
    fi <- apply(f,1,sum)
    fj <- apply(f,2,sum)
    ##    Dr = diag(fi)
    ##    Dc = diag(fj)
    f  <- (1/fi) * t(t(f)/fj)
    acp(f,wI=fi,wV=fj,center=TRUE,reduce=FALSE)
  }
