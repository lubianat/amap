###################################################
### chunk number 1: 
###################################################
cat("\n------  CLUSTERING TOOLS -------\n")
cat("\n------  Hierarchical clustering  -------\n")

data(USArrests)
h = hcluster(USArrests)
plot(h)
readline("Next")


###################################################
### chunk number 2: 
###################################################
cat("\n------  Hierarchical clustering using function heatmap -------\n")
heatmap(as.matrix(USArrests),
        hclustfun=hcluster,
        distfun=function(u){u})
readline("Next")


###################################################
### chunk number 3: 
###################################################
cat("\n------  Parralelized Hierarchical clustering  -------\n")
h = hclusterpar(USArrests,nbproc=4)   
readline("Next")


###################################################
### chunk number 4: 
###################################################
cat("\n------  K-means clustering  -------\n")
Kmeans(USArrests,centers=3,method="correlation")
readline("Next")


###################################################
### chunk number 5: 
###################################################
cat("\n------  ROBUST TOOLS  -------\n")
cat("\n------  A robust variance computation  -------\n")
data(lubisch)
lubisch <- lubisch[,-c(1,8)]
varrob(scale(lubisch),h=1)
readline("Next")


###################################################
### chunk number 6: 
###################################################
cat("\n------  A robust principal component analysis  -------\n")
p <- acpgen(lubisch,h1=1,h2=1/sqrt(2))
plot(p)
readline("Next")

###################################################
### chunk number 6: 
###################################################
cat("\n------  Another robust principal component analysis  -------\n")
p <- acprob(lubisch,h=4)
plot(p)
readline("Next")


cat("\n------ BUILDING HIERARCHICAL CLUSTERING WITH ANOTHER SOFTWARE  -------\n")
cat("\n------ Write data table to Xcluster file format   -------\n")
r2xcluster(USArrests,file='USArrests_xcluster.txt')
readline("Next")

cat("\n------ Write data table to cluster file format   -------\n")
r2cluster(USArrests,file='USArrests_cluster.txt')
readline("Next")

cat("\n------ Hierarchical clustering (need Xcluster tool  by Gavin Sherlock)  -------\ntry:\n\nh.xcl=xcluster(USArrests)\nplot(h.xcl)
") 



readline("Next")



hr = hcluster(USArrests)
hc = hcluster(t(USArrests))
cat("\n------  USING OTHER VISUALIZATION SOFTWARES  -------\n")

cat("\n------  Export hclust objects to Newick format files  -------\n")
write(hc2Newick(hr),file='hclust.newick')
readline("Next")

cat("\n------  Export hclust objects to Freeview or Treeview
visualization softwares  -------\n")
r2atr(hc,file="cluster.atr")
r2gtr(hr,file="cluster.gtr")
r2cdt(hr,hc,USArrests ,file="cluster.cdt")
readline("Next")

cat("\n------  Clustering and Export hclust objects to Freeview or Treeview
visualization softwares  -------\n")

hclust2treeview(USArrests,file="cluster.cdt")

