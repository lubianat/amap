r2xcluster <- function(data,labels=FALSE,description=FALSE,file="xcluster.txt"){
#-------------------------------------------------------
#
#  Created       : 10/12/01
#  Last Modified : Time-stamp: <2002-11-12 13:17:15 lucas>
#
#  Description   : Write to Xcluster file format (Xcluster make
#                  Hierarchical cluster analysis)
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#  See Also      : xcluster2r
#
#  Xcluster      : Xcluster is a C program that performs
#                  hierarchical clustering, K-means and
#                  SOM. 
#                  Xcluster is copyrighted. 
#                  To get or have information on Xcluster:
#          
#                  http://genome-www.stanford.edu/~sherlock/cluster.html     
#
#  Licence       : GPL (r2xcluster, not Xcluster)
#
#-------------------------------------------------------
#
#
# Example:
# source('r2xcluster.R')
# r2xcluster(data)
# system('Xcluster -f xcluster.txt')
# h <- xcluster2r('xcluster.gtr')
# library(mva)
# plot(h,hang=-1)
#
#-------------------------------------------------------


data <- as.data.frame(cbind(1,data))

# If the name column does not exist -> creation
# of this column else: put column name on
# first place.
if(!description){
        if(!labels){
                data <- cbind(row.names(data),data)
        }
        else{
                data <- cbind(as.factor(data[,2]),data)
        }
}
else{
	data <- cbind(as.factor(data[,3]),data[,-3])        
}
if(!labels){
        data <- cbind(row.names(data),data)
}
else{
        data <- cbind(as.factor(data[,3]),data[,-3])
}

levels(data[,1]) <- c(levels(data[,1]),"1","NAME","EWEIGHT")
levels(data[,2]) <- c(levels(data[,2]),"1","DESCRIPTION")
data <- rbind(1,data)
nom <- c("NAME","DESCRIPTION","GWEIGHT",names(data)[-c(1,2,3)])
data <- rbind(nom,data)
data[2,1] <- "EWEIGHT"

# NA will be replaced by "" in the file
#data[2:length(data[,1]),2] <- NA
data[2,2:3] <- NA

# Write file
write.table(data,file,sep='\t',row.names = FALSE,col.names = FALSE,na="",quote=FALSE)
}
