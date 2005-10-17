r2cluster <- function(data,labels=FALSE,colname="ACC",description=FALSE,file="cluster.txt",dec='.'){
#-------------------------------------------------------
#
#  Created       : 05/07/02
#  Last Modified : Time-stamp: <2003-07-22 16:53:00 lucas>
#
#  Description   : Write to Cluster file format (Cluster make
#                  Hierarchical cluster analysis)
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#  See Also      : r2xcluster
#
#  Xcluster      : Cluster is a program that performs
#                  hierarchical clustering, K-means and
#                  SOM. 
#                  Cluster is copyrighted. 
#                  To get or have information on Cluster:
#          
#                  http://rana.lbl.gov/EisenSoftware.htm
#
#  Licence       : GPL 
#
#
#-------------------------------------------------------
#
#
# Example:
# source('r2cluster.R')
# r2cluster(data)
#
#-------------------------------------------------------

n <- length(data[,1])
data <- as.data.frame(cbind(1,1:n,data))

# If the name column does not exist -> creation
# of this column, else: put column name on
# first place.
if(!description){
        if(!labels){
                data <- cbind(row.names(data),data)
        }
        else{
                data <- cbind(as.factor(data[,3]),data)
        }
}
else{
        data <- cbind(as.factor(data[,4]),data[,-4])
}
if(!labels){
        data <- cbind(row.names(data),data)
}
else{
        data <- cbind(as.factor(data[,4]),data[,-4])
}
levels(data[,1]) <- c(levels(data[,1]),"1",colname,"EWEIGHT")
levels(data[,2]) <- c(levels(data[,2]),"1","NAME","EWEIGHT")
data <- rbind(1,data)
nom <- c(colname,"NAME","GWEIGHT","GORDER",names(data)[-c(1,2,3,4)])
names(data) <- nom
data[1,1] <- "EWEIGHT"
data[1,2:4] <- NA

# NA will be replaced by "" in the file
#data[2:length(data[,1]),2] <- NA
#data[2,3] <- NA

# Write file
if(dec==',')
  {
    data<-apply(data,2,function(u){chartr(".",",",u)})
    data[data=="NA"]<-NA
  }

write.table(data,file,sep='\t',row.names = FALSE,col.names = TRUE,na="",quote=FALSE)
}
