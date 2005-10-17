xcluster <- function(data,distance="euclidean",clean=FALSE,tmp.in="tmp.txt",tmp.out="tmp.gtr"){
#-------------------------------------------------------
#
#  Created       : 10/12/01
#  Last Modified : Time-stamp: <2002-11-12 09:48:10 lucas>
#
#  Description   : Execute Xcluster (Xcluster make
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
#
#-------------------------------------------------------
#  Xcluster flags
#-------------------------------------------------------
#
# -g 0|1|2 ; 0 indicates no gene clustering, 1 indicates non-centered metric, 2 indicates centered metric when clustering genes.  2 is default.
# -e 0|1|2 ; 0 indicates no experiment clustering, see above for 1 and 2.  0 is the default.
# -p 0|1 ; whether to use pearson correlation (1), or Euclidean distance (0).  1 is the default.
# -s 0|1 ; whether to make a SOM (1). 0 is the default. 
# -x specify x dimension of SOM
# -y specify y dimension of SOM
# -r 0|1 ; whether to seed the random number generator with the time when making a SOM.  1 is the default.
# -k num ; how many k-means clusters to make.  num is an integer, greater than 1, and preferably less than a gazillion. 
# -l 0|1 ; whether to log transform the data.
# -u string ; a unique identifier by which to name the output files, instead of basing their names on the input file.  eg -u 888 will generate 888.cdt as an output file.
#-------------------------------------------------------
#
#
# Example:
# source('r2xcluster.R')
# system('Xcluster -f xcluster.txt')
# h <- xcluster2r('xcluster.gtr')
# library(mva)
# plot(h,hang=-1)
#
#-------------------------------------------------------

r2xcluster(data,file=tmp.in)

#-------------------------------------
#      CASE DISTANCE=EUCLIDEAN
#-------------------------------------
if(substr(distance,1,1)=="e")
   {
    script <- paste ("Xcluster -f",tmp.in," -e 0 -p 0 -s 0 -l 0")
    system(script)
    tree <- xcluster2r(file=tmp.out,distance="euclidean",labels=TRUE,fast= TRUE,clean=clean)
    }
#-------------------------------------
#      CASE DISTANCE=PEARSON
#-------------------------------------
if(substr(distance,1,1)=="p")
   {
    script <- paste ("Xcluster -f",tmp.in,"-g 2 -e 0 -p 1 -s 0 -l 0")
    system(script)
    tree <- xcluster2r(file=tmp.out,distance="pearson",labels=TRUE,fast= TRUE,clean=clean)
    }
#-------------------------------------
#   CASE DISTANCE=NONCENTEREDPEARSON
#-------------------------------------
if(substr(distance,1,1)=="n")
   {
    script <- paste ("Xcluster -f",tmp.in,"-g 1 -e 0 -p 1 -s 0 -l 0")
    system(script)
    tree <- xcluster2r(file=tmp.out,distance="pearson",labels=TRUE,fast= TRUE,clean=clean)
    }


script <-  paste ("rm ",substr(tmp.in,0,nchar(tmp.in)-3),"*",sep='')
system(script)
tree
}
