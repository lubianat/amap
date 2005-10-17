r2gtr <- function(hr,file="cluster.gtr",distance=hr$dist.method,dec=".",
                  digits=5)
{
#-------------------------------------------------------
#
#  Created       : 20/11/02	
#  Last Modified : Time-stamp: <2005-09-27 11:30:06 lucas>
#
#  Description   : Write hclust object to gtr atr
#                  files (Xcluster or Cluster output).
#                  Visualisation of cluster can be
#                  done with tools like treeview
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#  See Also      : plot.hclust  (library mva)
#
#  Licence       : GPL 
#-------------------------------------------------------



  height <- hr$height

  if(substr(distance,1,1)=="p")
    {
      height  <- 1 - height
    }
  else 
    {
      height <- height +1
      height  <- height[1] / height
#      height  <- (( 2 / height ) -1 ) / (( 2 / height[1] ) -1)
    }
  height  <- signif(height, digits = digits)
  
  
  n <- length(height)
  node <- 1:n
  node <- paste ('NODE',node,'X',sep='')

  merge1  <- hr$merge[,1]
  merge11 <- paste ('NODE',merge1,'X',sep='')
  merge12 <- paste ('GENE',-1-merge1,'X',sep='')
  merge1[hr$merge[,1]>0] <- merge11[hr$merge[,1]>0]
  merge1[hr$merge[,1]<0] <- merge12[hr$merge[,1]<0]
  
  merge2  <- hr$merge[,2]
  merge11 <- paste ('NODE',merge2,'X',sep='')
  merge12 <- paste ('GENE',-1-merge2,'X',sep='')
  merge2[hr$merge[,2]>0] <- merge11[hr$merge[,2]>0]
  merge2[hr$merge[,2]<0] <- merge12[hr$merge[,2]<0]


  data  <- data.frame(cbind(node,merge1,merge2))
  data  <- cbind(data,height)
 
  write.table(data,file=file,row.name=FALSE,col.names=FALSE,quote=FALSE,sep='\t',dec=dec)

}


#-----------------------------------
# Cosmetic modifications for r2atr
#-----------------------------------
r2atr <- function(hc,file="cluster.atr",distance=hc$dist.method,dec=".",
                  digits=5)
{



  height <- hc$height

  if(substr(distance,1,1)=="p")
    {
      height  <- 1 - height
    }
  else 
    {
      height <- height +1
      height  <- height[1] / height
    }
  height  <- signif(height, digits = digits)
  
  
  n <- length(height)
  node <- 1:n
  node <- paste ('NODE',node,'X',sep='')

  merge1  <- hc$merge[,1]
  merge11 <- paste ('NODE',merge1,'X',sep='')
  merge12 <- paste ('ARRY',-1-merge1,'X',sep='')
  merge1[hc$merge[,1]>0] <- merge11[hc$merge[,1]>0]
  merge1[hc$merge[,1]<0] <- merge12[hc$merge[,1]<0]
  
  merge2  <- hc$merge[,2]
  merge11 <- paste ('NODE',merge2,'X',sep='')
  merge12 <- paste ('ARRY',-1-merge2,'X',sep='')
  merge2[hc$merge[,2]>0] <- merge11[hc$merge[,2]>0]
  merge2[hc$merge[,2]<0] <- merge12[hc$merge[,2]<0]


  data  <- data.frame(cbind(node,merge1,merge2))
  data  <- cbind(data,height)
 
  if(dec==',')
    {
      data<-apply(data,2,function(u){chartr(".",",",u)})
      data[data=="NA"]<-NA
    }

  write.table(data,file=file,row.name=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

}

