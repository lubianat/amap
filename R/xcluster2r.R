xcluster2r <- function(file,distance="euclidean",labels=FALSE,fast=FALSE,clean=FALSE,dec="."){
#-------------------------------------------------------
#
#  Created       : 26/11/01
#  Last Modified : Time-stamp: <2002-11-22 11:48:29 lucas>
#
#  Description   : Read Xcluster output (Xcluster make
#                  Hierarchical cluster analysis) to
#                  analyse and plot the result. 
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#  See Also      : plot.hclust  (library mva)
#
#  Xcluster      : Xcluster is a C program that performs
#                  hierarchical clustering, K-means and
#                  SOM. 
#                  Xcluster is copyrighted. 
#                  To get or have information on Xcluster:
#          
#                  http://genome-www.stanford.edu/~sherlock/cluster.html     
#
#  Licence       : GPL (xcluster2r, not Xcluster)
#
#
#-------------------------------------------------------
#
#
# Example:
# source('xcluster2r.R')
# h <- xcluster2r('data.gtr')
# library(mva)
# plot(h,hang=-1)
#
#
#------------------------------
#  Determin if it is a .gtr or
#   .atr file
#------------------------------
ext <- substr(file,nchar(file)-2,nchar(file)-2)
if( (ext=='a') || (ext=='A') )
{
        premierelettre  <- 'A'
        atr             <- TRUE
}
else
{
        premierelettre  <- 'G'
        atr             <- FALSE
}

#------------------------------
#         reading file
#------------------------------



if(atr){
        premierelettre  <- 'A'
       }
else   {
        premierelettre  <- 'G'
       }

data    <- read.table(file,sep='\t',dec=dec)
data1   <- as.character(data[,2])
i1      <- (substr(data1,1,1)==premierelettre)*(-1)
#si i1 -1 => Gene, 
#Si i1  0 => Cluster
data1   <-  as.integer(substr(data1,5,nchar(data1)-1))
data1   <-  data1*i1+i1+data1*(1+i1)
data2   <- as.character(data[,3])
i2      <- (substr(data2,1,1)==premierelettre)*(-1)
data2   <-  as.integer(substr(data2,5,nchar(data2)-1))
data2   <-  data2*i2+i2+data2*(1+i2)
if(substr(distance,1,1)=="e")
   {
   Hheight <- 2 / (data[,4] +1 )   #Distance = Euclidean 
   }
   else{
   Hheight <- 1-data[,4]           #Distance = pearson (centered or not)
   }
rm(data)


#----------------------------------------------------
#         Giving ordered distances
#----------------------------------------------------
if(clean){
    Hheight <- cummax(Hheight)
    }



#----------------------------------------------------
#         Reorganizing data (Not essential, 
#         but gives the same look as R)
#         Fast=T skip this part
#----------------------------------------------------
if(!fast){
for(i in 1:length(data1)){
        if(abs(data1[i])>abs(data2[i])){
                tmp      <- data2[i]
                data2[i] <- data1[i]
                data1[i] <- tmp}
        if(data1[i]>0 && data2[i]<0)
        {
                tmp2     <- data2[i]
                data2[i] <- data1[i]
                data1[i] <- tmp2
        }
}
}
#----------------------------------------------------
#          Order the clusters 
#----------------------------------------------------

liste <- list("1"=c(data1[1],data2[1]))
for(i in 2:length(data1)){
        if(data1[i]<0 && data2[i]<0){
                liste[[as.character(i)]] <- c(data1[i],data2[i])
        }
        if(data1[i]>0 && data2[i]<0)
        {
                tmp <- as.character(data1[i])
                liste[[as.character(i)]] <- c(data2[i],liste[[tmp]])
                liste[[tmp]] <- c()
         }
        if(data1[i]<0 && data2[i]>0)
        {
                tmp <- as.character(data2[i])
                liste[[as.character(i)]] <- c(data1[i],liste[[tmp]])
                liste[[tmp]] <- c()
        }
        if(data1[i]>0 && data2[i]>0)
        {
                tmp1 <- as.character(data1[i])
                tmp2 <- as.character(data2[i])
                if(data1[i]<=data2[i])
                {
                        liste[[as.character(i)]] <- c(liste[[tmp1]],liste[[tmp2]])
                }
                else
                {
                        liste[[as.character(i)]] <- c(liste[[tmp2]],liste[[tmp1]])
                }                
                liste[[tmp1]] <- c()
                liste[[tmp2]] <- c()
        }
}


#--------------------------------
#       Giving the outputs
#-------------------------------- 

Hmerge  <- cbind(as.integer(data1),as.integer(data2))
Horder  <- -liste[[as.character(length(data2))]]
Hlabels <- as.character(1:(length(data2)+1))
Hmethod <- "average"
#Hdist.method <- "Pearson"
#Hcall  <- "Xcluster"

#--------------------------------
#  Getting labels on .cdt file
#--------------------------------
if(labels){
	if(atr)
	{
                filetxt <- paste(substr(file,0,nchar(file)-4),'.txt',sep='')
                data    <- read.table(filetxt,nrows=1,sep='\t',dec=dec)
                Hlabels <- as.character(t(data[4:length(data)]))
	}
	else
	{
                filecdt <- paste(substr(file,0,nchar(file)-4),'.cdt',sep='')
                data    <- read.table(filecdt,skip=1,sep='\t',dec=dec)
                data1   <- as.character(data[,1])
                data1   <- as.integer(substr(data1,5,nchar(data1)-1))
                data2   <- as.character(data[,2])
                Hlabels[data1+1] <- data2
        }
}

tree    <- list(merge=Hmerge,height=Hheight,order=Horder,labels=Hlabels,method=Hmethod)
class(tree) <- "hclust"
tree
}




