###### Curry et al. (in review), Artificial Intelligence in Geoscience
######R Code for Clustering and PCA 
######Clinopyroxene

##################################################################
##Step 1: Load necessary packages

# install.packages("compositions") ## do not need to install after first time
library(compositions)
# install.packages("rstudioapi")   ## do not need to install after first time
library(rstudioapi)

##################################################################
##Step 2: Read in data (Example data input provided here are plagioclase data from Table S1)

#Set working directory
rm(list = ls())
setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path)))

## read in data
inputdata <- read.csv("CentralSanJuanMinChemCompilation_Clinopyroxene.csv")

##Choose which dataset you want to cluster. For this paper, I chose to cluster only data from the publicaions Curry et al. (2021a) and Curry et al. (2021b)
inputdata <- inputdata[c(which(inputdata$DataSource == "Curryetal2021a")),] 

###Choose specific elements to use in clustering 
data1 <- inputdata[,c("SiO2_wt","TiO2_wt","Al2O3_wt","FeO_wt","MgO_wt","MnO_wt","CaO_wt","Eruption")]  ## makes sure the names exactly math those in your column headers

##################################################################
##Step 3: Data transformation and normalization

### isometric log ratio transformation (ilr)
### choose the elements you want to transform with the ilr method
elements<- c("SiO2_wt","TiO2_wt","Al2O3_wt","FeO_wt","MgO_wt","MnO_wt","CaO_wt") ## makes sure the names exactly math those in your column headers

dataset.compositional <- acomp(data1[,elements])
dataset.ilr <- ilr(dataset.compositional)
data1_ilr<- cbind(data1, dataset.ilr)

#### Data Normalization using median and standard deviation
##make a function
normalise <- function(x){
  num2 <- x-median(x)
  denom2 <- sd(x)
  return(num2/denom2)
}
##normalize
data1_norm<- as.data.frame(lapply(data1_ilr[,c("V1","V2","V3","V4","V5","V6")], normalise)) ##the number of vectors here is 1 less than the number of elements chosen above
tran.norm<- data1_norm   ##transformed and normalized dataset

data1_fin<- cbind(data1, tran.norm) #add your normalized and transformed vectors to your geochemical dataset

##NOTE: here you can export data for extremely randomized trees
write.csv(x = data1_fin, file = "cpx_medianSD_4clust_forRandTrees.csv")

##################################################################
##Step 4: Data Clustering
dist_mat<- dist(tran.norm, method = "euclidean")
clusters_1<- hclust(dist_mat, method = "ward.D2") 

##Step 4a Choose number of clusters
###Note, the number of clusters is determined by: 1) PCA (see Step 6 below) and 2) plotting various binary geochemical plots to see how well the data differentiate into clusters (not included in this code).
best_n_clusters<- 4

##Step 4b: Data Organization
clusterCut<-matrix(ncol=10, nrow = nrow(tran.norm))
for (i in 1:10){
  clusterCut[,i]<- cutree(clusters_1, k = i)
}
data1_clusters <- data.frame(data1_fin[,])
data1_clusters$Cluster <- clusterCut[,best_n_clusters]  ##add clusters to dataframe

##NOTE: here you can export data with clusters
write.csv(x = data1_clusters, file = "cpx_medianSD_4clust_forRandTrees_wClusters.csv")

cluster1 <- which(data1_clusters$Cluster == 1)
cluster2 <- which(data1_clusters$Cluster == 2)
cluster3 <- which(data1_clusters$Cluster == 3)
cluster4 <- which(data1_clusters$Cluster == 4)

##Step 4c Cluster Dendrograms 
seq.clusters<- data1_clusters$Cluster[clusters_1$order]
seq.clusters1<-NULL
for (i in 1:length(seq.clusters)){
  seq.clusters1[i]<- seq.clusters[i]-seq.clusters[i+1]
}
seq.clusters2<- c(seq.clusters[which(seq.clusters1>0 | seq.clusters1<0)], seq.clusters[length(seq.clusters)-1])
dists<- clusters_1$height[(length(clusters_1$height)-(best_n_clusters-2)):length(clusters_1$height)]

par(mfrow=c(1,1))
plot(clusters_1)
rect.hclust(clusters_1 , k = best_n_clusters, border = seq.clusters2)
legend("topleft", legend=seq.clusters2, col=seq.clusters2, lty=1,cex=0.5)

##################################################################
##Step 5: Pie Charts#########

##all data, all eruptions
Cl1_perc <- (length(cluster1)/length(data1_clusters$Cluster))*100
Cl2_perc <- (length(cluster2)/length(data1_clusters$Cluster))*100
Cl3_perc <- (length(cluster3)/length(data1_clusters$Cluster))*100
Cl4_perc <- (length(cluster4)/length(data1_clusters$Cluster))*100
par(mfrow=c(1,1))
pie(c(Cl1_perc,Cl2_perc,Cl3_perc,Cl4_perc),labels = c("C1 - 26%","C2 - 59%","C3 - 6%", "C4 - 9%"), 
    col=c("black","red","green","blue"),main="Cpx Clusters - All")

##Snowshoe Mountain Tuff only
data1_clustersSM <- data1_clusters[which(data1_clusters$Eruption == "SnowshoeMountainTuff"),]
SM1 <- length(which(data1_clustersSM$Cluster == 1))/length(data1_clustersSM$Cluster)*100
SM2 <-length(which(data1_clustersSM$Cluster == 2))/length(data1_clustersSM$Cluster)*100
SM3 <-length(which(data1_clustersSM$Cluster == 3))/length(data1_clustersSM$Cluster)*100
SM4 <-length(which(data1_clustersSM$Cluster == 4))/length(data1_clustersSM$Cluster)*100
pie(c(SM1,SM2,SM3,SM4),labels = c("SM1-31","SM2-69","SM3-0", "C4 - 0%"), 
    col=c("black","red","green","blue"),main="Cpx Clusters - SM")

##Nelson Mountain Tuff only
data1_clustersNM <- data1_clusters[which(data1_clusters$Eruption == "NelsonMountainTuff"),]
NM1 <-length(which(data1_clustersNM$Cluster == 1))/length(data1_clustersNM$Cluster)*100
NM2 <-length(which(data1_clustersNM$Cluster == 2))/length(data1_clustersNM$Cluster)*100
NM3 <-length(which(data1_clustersNM$Cluster == 3))/length(data1_clustersNM$Cluster)*100
NM4 <-length(which(data1_clustersNM$Cluster == 4))/length(data1_clustersNM$Cluster)*100
pie(c(NM1,NM2,NM3,NM4),labels = c("NM1-39","NM2-61","NM3-0", "C4 - 0%"), 
    col=c("black","red","green","blue"),main="Cpx Clusters - NM")


##Cebolla Creek Tuff only
data1_clustersCC <- data1_clusters[which(data1_clusters$Eruption == "CebollaCreekTuff"),]
CC1 <-length(which(data1_clustersCC$Cluster == 1))/length(data1_clustersCC$Cluster)*100
CC2 <-length(which(data1_clustersCC$Cluster == 2))/length(data1_clustersCC$Cluster)*100
CC3 <-length(which(data1_clustersCC$Cluster == 3))/length(data1_clustersCC$Cluster)*100
CC4 <-length(which(data1_clustersCC$Cluster == 4))/length(data1_clustersCC$Cluster)*100
pie(c(CC1,CC2,CC3,CC4),labels = c("CC1-26","CC2-74","CC3-0", "C4 - 0%"), 
    col=c("black","red","green","blue"),main="Cpx Clusters - CC")


##Rat Creek Tuff only
data1_clustersRC <- data1_clusters[which(data1_clusters$Eruption == "RatCreekTuff"),]
RC1 <-length(which(data1_clustersRC$Cluster == 1))/length(data1_clustersRC$Cluster)*100
RC2 <-length(which(data1_clustersRC$Cluster == 2))/length(data1_clustersRC$Cluster)*100
RC3 <-length(which(data1_clustersRC$Cluster == 3))/length(data1_clustersRC$Cluster)*100
RC4 <-length(which(data1_clustersRC$Cluster == 4))/length(data1_clustersRC$Cluster)*100
pie(c(RC1,RC2,RC3,RC4),labels = c("RC1-0","RC2-14","RC3-35", "C4 - 51%"), 
    col=c("black","red","green","blue"),main="Cpx Clusters - RC")

##################################################################
##Step 6: PCA
data.pca <- prcomp(tran.norm, center = TRUE,scale. = TRUE)
summary.pca<- summary(data.pca)
summary.pca #This shows the comulative proportion of variance explained by the PCA componenets

PCA1<- data.pca$x[,1]
PCA2<- data.pca$x[,2]
data1_clusters$color1 <- "black"
data1_clusters$color1[cluster1] <- "black"
data1_clusters$color1[cluster2] <- "red"
data1_clusters$color1[cluster3] <- "green"
data1_clusters$color1[cluster4] <- "blue"
plot(PCA1, PCA2, col= data1_clusters$color1, pch=19, cex = 1.5)

