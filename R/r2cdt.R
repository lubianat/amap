r2cdt <- function(hr,hc,data,labels=FALSE,description=FALSE,file="cluster.cdt",dec="."){
#-------------------------------------------------------
#
#  Created       : 20/11/02
#  Last Modified : Time-stamp: <2003-07-22 16:52:28 lucas>
#
#  Description   : Write data object to cdt 
#                  file (Xcluster or Cluster output).
#                  Should be use with r2gtr
#                  Visualisation of cluster can be
#                  done with tools like treeview
#
#                  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#
#  Licence       : GPL (r2xcluster, not Xcluster)
#
#-------------------------------------------------------
#
#

  n <- dim(data)[1]

  # Add GWEIGHT column
  data <- as.data.frame(cbind(1,data))
  

  # If the name column does not exist -> creation
  # of this column else: put column name on
  # first place. (column NAME in the file)
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
  # column YORF in the file
  if(!labels){
    data <- cbind(row.names(data),data)
  }
  else{
    data <- cbind(as.factor(data[,3]),data[,-3])
  }

  # add GID column
  GID <- paste ('GENE',0:(n-1),'X',sep='')
  data <- cbind(as.factor(GID),data)

  # Put data in right order
  data <- data [hr$order,]

  m <- dim(data)[2]
  data[,5:m] <- data[,5:m][hc$order]
  colnames(data)[5:m] <-  colnames(data)[5:m][hc$order]
  
  # Round data
  data[,5:m] <- signif(data[,5:m], digits = 4)


  levels(data[,1]) <- c(levels(data[,1]),"1","GID","AID","EWEIGHT")
  levels(data[,2]) <- c(levels(data[,2]),"1","UNIQID")
  levels(data[,3]) <- c(levels(data[,3]),"1","NAME")
  data <- rbind(1,data)
  nom <- c("GID","UNIQID","NAME","GWEIGHT",names(data)[-c(1,2,3,4)])
  names(data) <- nom
  data[1,1] <- "EWEIGHT"
  data <- rbind(c("AID",NA,NA,NA,paste('ARRY',hc$order-1,'X',sep='')),data)
  
  # NA will be replaced by "" in the file
  # data[2:length(data[,1]),2] <- NA
  data[2,2:4] <- NA

  # Write file
  if(dec==',')
    {
      data<-apply(data,2,function(u){chartr(".",",",u)})
      data[data=="NA"]<-NA
    }
  write.table(data,file = file,sep = '\t',row.names = FALSE,col.names = TRUE,na="",quote=FALSE)

}
