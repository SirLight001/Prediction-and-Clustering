#######################################################################################################
# PROJECT LECTURER: PROF. ADAM ZAGDANSKI                                                              #
# COURSE TITLE: DATA MINING                                                                           #
# PROJECT TOPIC: CLUSTER ANALYSIS OF BREAST CANCER CELL: AN APPLICATION OF PARTITION AND              #
#               HIERARCHICAL APPROACH                                                                 #
# STUDENTS AND ID: SEGUN LIGHT JEGEDE (257389) and ISAAC AKOJI PAUL (257388)                          #
#######################################################################################################
#-------------------------------
## LOADING THE DATA INTO R
#----------------------------------------------------------------------------------------------------------------------------
#Reading the data into r and renaming the variable names in a readable manner.
col.names=c("id_number","clump_thickness","uniformity_cell_size","uniformity_cell_shape",
            "marginal_adhesion","single_epithelial_cell_size","bare_nuclei",
            "bland_chromatin","normal_nucleoli","mitoses","class")
bcw <- read.csv("C:/Users/jeged/Downloads/breast-cancer-wisconsin.data", header=FALSE, col.names=col.names)
#View(bcw)
attach(bcw)
#-------------------------------
## DATA PREPARATION AND CLEANING
#----------------------------------------------------------------------------------------------------------------------------
library(DataExplorer)

#Checking the data type of each columns
str(bcw)
bcw$class = as.factor(bcw$class) #convert the class to factor with "2" as benign and "4" as malignant
levels(bcw$class)[levels(bcw$class)=="2"] <- "benign"
levels(bcw$class)[levels(bcw$class)=="4"] <- "malignant"
bcw[,2:10] <- suppressWarnings(apply(bcw[, 2:10], 2, function(x) as.numeric(as.character(x)))) #format all features as numeric
bcw$id_number = as.character(bcw$id_number) #id_number is nothing but a string of cells identification number
str(bcw) #Every attributes is now in their respective perfect form.

#Handling missing data by deleting the corresponding rows if the missing observations are not too much
t(introduce(bcw))
sum(is.na(bcw)) #check for missing observations
plot_intro(bcw)
plot_missing(bcw) 
bcw<-na.omit(bcw) #2.29% of the bare_nuclei measurement are missing variables, thus, we decided to remove any form of missing observation.
sum(is.na(bcw))
plot_missing(bcw) 
nrow(bcw) #The data reduced from 699 to 683, we suppose we did not lose too much information, just about 3%
#View(bcw)

#-------------------------------
##EXPLORATORY DATA ANALYSIS
#----------------------------------------------------------------------------------------------------------------------------
#Describing the Grouping Variable
library(tidyverse)
library(scales)
bcwnew <- bcw %>% group_by(class) %>% 
  summarize(count = n()) %>%  # count records by species
  mutate(percentage = count/sum(count))  # find percent of total

ggplot(bcwnew, aes(class, percentage, fill = class)) + 
  geom_bar(stat='identity') + 
  geom_text(aes(label=scales::percent(percentage)), position = position_stack(vjust = .5))+
  scale_y_continuous(labels = scales::percent)

#Describing the features
my.summary <- function(df)
{
  results <- matrix(, nrow = 9, ncol = ncol(df))
  for (i in 1:ncol(df)){
    X=df[,i]
    results[,i] <- rbind(min(X),quantile(X,0.25), median(X), mean(X), quantile(X,0.75), max(X), var(X), sd(X), IQR(X))
  }
  rownames(results) <- c("min", "Q1", "median", "mean", "Q3", "max", "var", "sd", "IQR")
  colnames(results) <-names(df)
  return(results)
}
ms<-my.summary(bcw[,2:10]) 
ms
write.table(ms, file = "summary statistics.txt", sep = ",", quote = FALSE, row.names = F)

#construct the plots three by three
#construct the histogram plots
library(ggpubr)
ha<-fg<-ggplot(bcw, aes(x = clump_thickness, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hb<-ggplot(bcw, aes(x = uniformity_cell_size, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hc<-ggplot(bcw, aes(x = uniformity_cell_shape, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hd<-fg<-ggplot(bcw, aes(x = marginal_adhesion, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
ggarrange(ha,hb,hc,hd,labels = c("A", "B", "C","D"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

hd<-fg<-ggplot(bcw, aes(x = marginal_adhesion, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
he<-ggplot(bcw, aes(x = single_epithelial_cell_size, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hf<-ggplot(bcw, aes(x = bare_nuclei, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hg<-fg<-ggplot(bcw, aes(x = bland_chromatin, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
hh<-ggplot(bcw, aes(x = normal_nucleoli, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
ggarrange(he,hf,hg,hh,labels = c("E", "F","G", "H"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

hi<-ggplot(bcw, aes(x = mitoses, fill = class)) + geom_histogram(position = "identity", alpha = 0.4)
ggarrange(hi,labels = c("I"),ncol = 2, nrow = 2) #, common.legend = TRUE, legend="bottom")

'plot_histogram(bcw)'

#construct the density plots
plot_density(bcw[2:4])
plot_density(bcw[5:7])
plot_density(bcw[8:10])

#construct the normal qq plot
plot_qq(bcw[2:4])
plot_qq(bcw[5:7])
plot_qq(bcw[8:10])

#construct the barplots
bpa <- ggplot(bcw, aes(x = class, y = clump_thickness))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpb <- ggplot(bcw, aes(x = class, y = uniformity_cell_size))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpc <- ggplot(bcw, aes(x = class, y = uniformity_cell_shape))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpd <- ggplot(bcw, aes(x = class, y = marginal_adhesion))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
ggarrange(bpa,bpb,bpc,bpd,labels = c("A", "B", "C","D"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

bpe <- ggplot(bcw, aes(x = class, y = single_epithelial_cell_size))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpf <- ggplot(bcw, aes(x = class, y = bare_nuclei))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bpg <- ggplot(bcw, aes(x = class, y = bland_chromatin))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
bph <- ggplot(bcw, aes(x = class, y = normal_nucleoli))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
ggarrange(bpe,bpf,bpg,bph,labels = c("E", "F", "G","H"),ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

bpi <- ggplot(bcw, aes(x = class, y = mitoses))+geom_boxplot(aes(color = class))+
  scale_color_manual(values = c("#00AFBB", "#E7B800"))
ggarrange(bpi,labels = c("I"),ncol = 2, nrow = 2) #, common.legend = TRUE, legend="bottom")

'plot_boxplot(bcw, by="class")'

#Construct Boxplot without grouping
library(reshape)
bcwData <- melt(bcw)
par(mar=c(10,7,1,1))
boxplot(data=bcwData, value~variable, las=2)
#pairs(bcw[2:10], pch = 21, bg = c("#d95f02", "#7570b3")[unclass(bcw$class)])

plot_correlation(bcw, type = "continuous") #correlation plot
#we can also do sum data preparation here; it is an error if the measurement does not fall within 1-10 i.e Minimum and Maximum

#-------------------------------
## NORMALIZATION
#----------------------------------------------------------------------------------------------------------------------------
a<-(bcw$clump_thickness-min(bcw$clump_thickness))/(max(bcw$clump_thickness)-min(bcw$clump_thickness))
b<-(bcw$uniformity_cell_size-min(bcw$uniformity_cell_size))/(max(bcw$uniformity_cell_size)-min(bcw$uniformity_cell_size))
c<-(bcw$uniformity_cell_shape-min(bcw$uniformity_cell_shape))/(max(bcw$uniformity_cell_shape)-min(bcw$uniformity_cell_shape))
d<-(bcw$marginal_adhesion-min(bcw$marginal_adhesion))/(max(bcw$marginal_adhesion)-min(bcw$marginal_adhesion))
e<-(bcw$single_epithelial_cell_size-min(bcw$single_epithelial_cell_size))/(max(bcw$single_epithelial_cell_size)-min(bcw$single_epithelial_cell_size))
f<-(bcw$bare_nuclei-min(bcw$bare_nuclei))/(max(bcw$bare_nuclei)-min(bcw$bare_nuclei))
g<-(bcw$bland_chromatin-min(bcw$bland_chromatin))/(max(bcw$bland_chromatin)-min(bcw$bland_chromatin))
h<-(bcw$normal_nucleoli-min(bcw$normal_nucleoli))/(max(bcw$normal_nucleoli)-min(bcw$normal_nucleoli))
i<-(bcw$mitoses-min(bcw$mitoses))/(max(bcw$mitoses)-min(bcw$mitoses))
bcw1<-data.frame(bcw[,1],a,b,c,d,e,f,g,h,i,bcw[,11])
colnames(bcw1)<-c("ID","clump_thickness","uniformity_cell_size","uniformity_cell_shape","marginal_adhesion",
                  "single_epithelial_cell_size","bare_nuclei","bland_chromatin","normal_nucleoli","mitoses","class")
str(bcw1)
#View(bcw1)

#Boxplot of the normalized data
bcw1Data <- melt(bcw1)
par(mar=c(10,7,1,1))
boxplot(data=bcw1Data, value~variable, las=2)


#-------------------------------
## Clustering Start
#----------------------------------------------------------------------------------------------------------------------------
library(stats)
library(cluster)
library("factoextra")
library(rgl)
library(scatterplot3d)
bcw.features <- bcw[,2:10] # We remove class labels
bcw.real.class.labels <- bcw[,11]
bcw.names <- paste(bcw$class, bcw$id_number, sep=" ") #we assign label names by including ID

#-------------------------------
## ACCESSING CLUSTERING TENDENCY
#----------------------------------------------------------------------------------------------------------------------------
# Random data generated from the iris data set
random_df <- apply(bcw.features, 2, function(x){runif(length(x), min(x), (max(x)))})
random_df <- as.data.frame(random_df)

fviz_pca_ind(prcomp(bcw.features), title = "PCA - Breast Cancer data", habillage = bcw.real.class.labels,  palette = "jco",
             geom = "point", ggtheme = theme_classic(),legend = "bottom") # Plot bcw data set

fviz_pca_ind(prcomp(random_df), title = "PCA - Random data", geom = "point", ggtheme = theme_classic()) # Plot the random df

res.bcw.real <- get_clust_tendency(bcw.features, n = nrow(bcw.features)-1, graph = FALSE)
res.bcw.real$hopkins_stat

res.bcw.random <- get_clust_tendency(random_df, n = nrow(random_df)-1, graph = FALSE)
res.bcw.random$hopkins_stat

#-----------------------------------------
## Selecting the Optimal Number of Cluster
#----------------------------------------------------------------------------------------------------------------------------
# Naive approach: 'elbow method' - We are looking for a strong bend in the chart, the so-called "elbow" or "knee".
fviz_nbclust(bcw.features, FUNcluster = kmeans, method="wss", k.max=10) + geom_vline(xintercept=2, linetype=2) #KMeans
fviz_nbclust(bcw.features, FUNcluster = cluster::pam, method="wss", k.max=10) + geom_vline(xintercept=2, linetype=2) #PAM
fviz_nbclust(bcw.features, FUNcluster = cluster::clara, method="wss", k.max=10) + geom_vline(xintercept=2, linetype=2) #CLARA
fviz_nbclust(bcw.features, FUNcluster = hcut, method="wss", k.max=10) + geom_vline(xintercept=2, linetype=2) # hierarchical clustering

# Other advanced methods used to select the optimal K: Silhouette.
fviz_nbclust(bcw.features, FUNcluster = kmeans, method = "silhouette") #KMeans
fviz_nbclust(bcw.features, FUNcluster = cluster::pam, method = "silhouette") # PAM
fviz_nbclust(bcw.features, FUNcluster = cluster::clara, method = "silhouette") # CLARA
fviz_nbclust(bcw.features, FUNcluster = hcut, method = "silhouette") # hierarchical clustering

# Using the NbClust
library(NbClust)
NbClust.results.1 <- NbClust(bcw.features, distance="euclidean", min.nc=2, max.nc=10, method="complete", index="all")
NbClust.results.1$All.index
NbClust.results.1$Best.nc
NbClust.results.1$Best.partition
factoextra::fviz_nbclust(NbClust.results.1) + theme_minimal() + ggtitle("Optimal number of clusters")

#-------------------------------
## Internal Cluster Validation
#----------------------------------------------------------------------------------------------------------------------------
library(clValid)
library(mclust)
'methods <- c("agnes","kmeans", "diana", "pam", "clara")
K.range <- 2:5 # range for number of clusters
internal.validation <- clValid(bcw.features, nClust=K.range, clMethods=methods, validation="internal")
y
summary(internal.validation)
optimalScores(internal.validation)
par(mfrow = c(2, 2))
plot(internal.validation, legend = FALSE, lwd=2)
plot.new()
legend("center", clusterMethods(internal.validation), col=1:9, lty=1:9, pch=paste(1:9))

stability.validation <- clValid(bcw.features, nClust=K.range, clMethods=methods, validation="stability")
y
summary(stability.validation)
optimalScores(stability.validation)
par(mfrow = c(2,2))
plot(stability.validation, measure=c("APN","AD","ADM"), legend=FALSE, lwd=2)
plot.new()
legend("center", clusterMethods(stability.validation), col=1:9, lty=1:9, pch=paste(1:9))'

#-------------------------------
## PCA
#----------------------------------------------------------------------------------------------------------------------------
bcw.pca <- bcw[, 2:10]
prcomp(bcw.pca, retx=T, center=T, scale.=T) -> bcw.after.pca

library(factoextra)
fviz_eig(bcw.after.pca) #The Knee Plot for PCA 
fviz_pca_var(bcw.after.pca, col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Results for Variables
res.bcw.pca <- get_pca_var(bcw.after.pca)
res.bcw.pca$coord          # Coordinates
res.bcw.pca$contrib        # Contributions to the PCs
res.bcw.pca$cos2           # Quality of representation

library("corrplot")
corrplot(res.bcw.pca$cos2, is.corr=FALSE)
fviz_contrib(bcw.after.pca, choice = "var", axes = 1:2, top = 9) # Contributions of variables to PC1 and PC2
fviz_contrib(bcw.after.pca, choice = "var", axes = 1, top = 9) # Contributions of variables to PC1
fviz_contrib(bcw.after.pca, choice = "var", axes = 2, top = 9) # Contributions of variables to PC2

print("Principal components:")
print(bcw.after.pca$rotation)
summary(bcw.after.pca) #summary of loadings on components

# We analyse the amount of variance explained by subsequent principal components
variance <- round(((bcw.after.pca$sdev ^2)/sum(bcw.after.pca$sdev^2)), 4)
cumulative.variance <- cumsum(variance)
pca.df.var <- data.frame(PCs=c("Dim_1","Dim_2","Dim_3","Dim_4","Dim_5","Dim_6","Dim_7","Dim_8","Dim_9"),
                         Variance=variance, Cummulative_Variance=cumulative.variance)
ggplot(pca.df.var, aes(x=PCs, y=Variance, fill=PCs)) + geom_bar(stat="identity") + 
  geom_text(aes(label=Variance), vjust=-0.3, size=3.5) #Variance Explained by PCA
ggplot(pca.df.var, aes(x=PCs, y=Cummulative_Variance, fill=PCs)) + geom_bar(stat="identity") + 
  geom_text(aes(label=Cummulative_Variance), vjust=-0.3, size=3.5) #Cummulative Variance by PCA

pca.features<-data.frame(bcw.after.pca$x[,1], bcw.after.pca$x[,2], bcw.after.pca$x[,3])
colnames(bcw1)<-c("PC1","PC2","PC3")

#-------------------------------
## K-MEANS
#----------------------------------------------------------------------------------------------------------------------------
k <- 2 # Partition into K clusters

### FOR ALL FEATURES
kmeans.k2.10x <- kmeans(bcw.features, 2, iter.max=10, nstart=10)
bcw.kmeans.labels <- kmeans.k2.10x$cluster
plot(bcw.features, col=kmeans.k2.10x$cluster)
title('K-Means Clustering for Breast Cancer Problem (10 random initialization)')

bcw.sil.kmeans <- silhouette(bcw.kmeans.labels, dist(bcw.features))
fviz_silhouette(bcw.sil.kmeans, xlab="K-means") #silhouette information

# Visualization of cluster analysis results (3D scatterplot)
plot3d(bcw.features$uniformity_cell_size, bcw.features$uniformity_cell_shape, bcw$clump_thickness, col=bcw.kmeans.labels, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="uniformity_cell_size", 
       ylab="uniformity_cell_shape", zlab="clump thickness")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = bcw.kmeans.labels, cex=1, inset=c(0.02))


### FOR SELECTED PCA FEATURES
kmeans.k2.10xpca <- kmeans(pca.features, centers=k, iter.max=10, nstart=10)
pca.kmeans.labels <- kmeans.k2.10xpca$cluster
plot(pca.features, col=kmeans.k2.10xpca$cluster)
title('K-Means Clustering for Breast Cancer Problem PCA (10 random initialization)')

pca.sil.kmeans <- silhouette(pca.kmeans.labels, dist(pca.features))
fviz_silhouette(pca.sil.kmeans, xlab="K-means") #silhouette information

plot3d(pca.features$bcw.after.pca.x...1., pca.features$bcw.after.pca.x...2., pca.features$bcw.after.pca.x...3., col=pca.kmeans.labels, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="Dim_1", 
       ylab="Dim_2", zlab="Dim_3")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = pca.kmeans.labels, cex=1, inset=c(0.02))

#-------------------------------
## PAM
#----------------------------------------------------------------------------------------------------------------------------
library(cluster)

### Application of PAM algorithm FOR ALL FEATURES
bcw.pam2 <- pam(x=bcw.features, k=2)
X11()
plot(bcw.pam2) # default visualization (note: plot() works differently for quantitative and mixed data types)
(summary(bcw.pam2)) 
bcw.cluster.labels <- bcw.pam2$clustering

bcw.sil.pam2 <- silhouette(bcw.cluster.labels, dist(bcw.features))
fviz_silhouette(bcw.sil.pam2, xlab="PAM") #silhouette information

# Visualization of cluster analysis results (3D scatterplot)
plot3d(bcw.features$uniformity_cell_size, bcw.features$uniformity_cell_shape, bcw$clump_thickness, col=bcw.cluster.labels, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="uniformity_cell_size", 
       ylab="uniformity_cell_shape", zlab="clump thickness")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = bcw.cluster.labels, cex=1, inset=c(0.02))
#snapshot3d(filename = '3dplot.png', fmt = 'png')


### Application of PAM algorithm FOR PCA FEATURES
pca.pam2 <- pam(x=pca.features, k=2)
#X11()
#plot(pca.pam2) # default visualization (note: plot() works differently for quantitative and mixed data types)
(summary(pca.pam2)) 
pca.cluster.labels <- pca.pam2$clustering

pca.sil.pam2 <- silhouette(pca.cluster.labels, dist(pca.features))
fviz_silhouette(pca.sil.pam2, xlab="PAM") #silhouette information

# Visualization of cluster analysis results (3D scatterplot)
plot3d(pca.features$bcw.after.pca.x...1., pca.features$bcw.after.pca.x...2., pca.features$bcw.after.pca.x...3., col=pca.cluster.labels, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="Dim_1", 
       ylab="Dim_2", zlab="Dim_3")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = pca.cluster.labels, cex=1, inset=c(0.02))
#snapshot3d(filename = '3dplot.png', fmt = 'png')


#-------------------------------
## CLARA
#----------------------------------------------------------------------------------------------------------------------------
# compute CLARA FOR ALL FEATURES
bcw.clara <- clara(bcw.features, 2, samples=200, pamLike = TRUE)
print(bcw.clara)
bcw.clara.clust <- bcw.clara$cluster

bcw.sil.clara <- silhouette(bcw.clara.clust, dist(bcw.features))
fviz_silhouette(bcw.sil.clara, xlab="CLARA")

# Visualization of cluster analysis results (3D scatterplot)
plot3d(bcw.features$uniformity_cell_size, bcw.features$uniformity_cell_shape, bcw$clump_thickness, col=bcw.clara.clust, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="uniformity_cell_size", 
       ylab="uniformity_cell_shape", zlab="clump thickness")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = bcw.clara.clust, cex=1, inset=c(0.02))
#snapshot3d(filename = '3dplot.png', fmt = 'png')


# compute CLARA FOR PCA FEATURES
pca.clara <- clara(pca.features, 2, samples=200, pamLike = TRUE)
print(pca.clara)
pca.clara.clust <- pca.clara$cluster

pca.sil.clara <- silhouette(pca.clara.clust, dist(pca.features))
fviz_silhouette(pca.sil.clara, xlab="CLARA")

# Visualization of cluster analysis results (3D scatterplot)
plot3d(pca.features$bcw.after.pca.x...1., pca.features$bcw.after.pca.x...2., pca.features$bcw.after.pca.x...3., col=pca.clara.clust, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="Dim_1", 
       ylab="Dim_2", zlab="Dim_3")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = pca.clara.clust, cex=1, inset=c(0.02))
#snapshot3d(filename = '3dplot.png', fmt = 'png')

#-------------------------------
## AGNES
#----------------------------------------------------------------------------------------------------------------------------

## FOR ALL FEATURES
# We compare the methods available for AGNES
m <- c( "average", "single", "complete")
names(m) <- c( "average", "single", "complete")
# function to compute coefficient
ac <- function(x) {
  agnes(bcw.features, method = x)$ac
}
map_dbl(m, ac)  

bcw.agnes.avg <- agnes(bcw.features, method = "average")
pltree(bcw.agnes.avg, cex = 0.6, hang = -1, main = "Dendrogram of AGNES with Average Linkage")
(bcw.agnes.avg.k2 <- cutree(bcw.agnes.avg, k=2)) #Cutting off at k=2
table(bcw.agnes.avg.k2)

bcw.sil.agnes <- silhouette(bcw.agnes.avg.k2, dist(bcw.features))
fviz_silhouette(bcw.sil.agnes, xlab="AGNES")

fviz_dend(bcw.agnes.avg, cex=0.4, main="Dendrogram of AGNES with Average Linkage") # standard dendrogram
fviz_dend(bcw.agnes.avg, k=2, cex=0.4) # clustered dendrogram
fviz_dend(bcw.agnes.avg, type="circular", cex=0.4, k=2,  main="Dendrogram of AGNES with Average Linkage") # circular dendrogram

# Visualization of cluster analysis results (3D scatterplot)
plot3d(bcw.features$uniformity_cell_size, bcw.features$uniformity_cell_shape, bcw$clump_thickness, col=bcw.agnes.avg.k2, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="uniformity_cell_size", 
       ylab="uniformity_cell_shape", zlab="clump thickness")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = bcw.agnes.avg.k2, cex=1, inset=c(0.02))


## FOR PCA FEATURES
# We compare the methods available for AGNES
m <- c( "average", "single", "complete")
names(m) <- c( "average", "single", "complete")
# function to compute coefficient
ac <- function(x) {
  agnes(pca.features, method = x)$ac
}
map_dbl(m, ac)  

pca.agnes.avg <- agnes(pca.features, method = "average")
pltree(pca.agnes.avg, cex = 0.6, hang = -1, main = "Dendrogram of AGNES with Average Linkage - PCA")
(pca.agnes.avg.k2 <- cutree(pca.agnes.avg, k=2)) #Cutting off at k=2
table(pca.agnes.avg.k2)

pca.sil.agnes <- silhouette(pca.agnes.avg.k2, dist(pca.features))
fviz_silhouette(pca.sil.agnes, xlab="AGNES")

fviz_dend(pca.agnes.avg, cex=0.4, main="Dendrogram of AGNES with Average Linkage - PCA") # standard dendrogram
fviz_dend(pca.agnes.avg, k=2, cex=0.4) # clustered dendrogram
fviz_dend(pca.agnes.avg, type="circular", cex=0.4, k=2,  main="Dendrogram of AGNES with Average Linkage - PCA") # circular dendrogram

# Visualization of cluster analysis results (3D scatterplot)
plot3d(pca.features$bcw.after.pca.x...1., pca.features$bcw.after.pca.x...2., pca.features$bcw.after.pca.x...3., col=pca.agnes.avg.k2, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="Dim_1", 
       ylab="Dim_2", zlab="Dim_3")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = pca.agnes.avg.k2, cex=1, inset=c(0.02))

#-------------------------------
## DIANA
#----------------------------------------------------------------------------------------------------------------------------

# compute divisive hierarchical clustering FOR ALL FEATURES
bcw.diana <- diana(bcw.features)
bcw.diana$dc
pltree(bcw.diana, cex = 0.6, hang = -1, main = "Dendrogram of diana")
rect.hclust(bcw.diana, k = 2, border = 2:10)
bcw.diana.clust <- cutree(bcw.diana, k = 2)

bcw.sil.diana <- silhouette(bcw.diana.clust, dist(bcw.features))
fviz_silhouette(bcw.sil.diana, xlab="DIANA")

fviz_cluster(list(data = bcw.features, cluster = bcw.diana.clust))
fviz_dend(bcw.diana, k=2, cex=0.4)
fviz_dend(bcw.diana, type="circular", cex=0.4, k=2)

# Visualization of cluster analysis results (3D scatterplot)
plot3d(bcw.features$uniformity_cell_size, bcw.features$uniformity_cell_shape, bcw$clump_thickness, col=bcw.diana.clust, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="uniformity_cell_size", 
       ylab="uniformity_cell_shape", zlab="clump thickness")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = bcw.diana.clust, cex=1, inset=c(0.02))


# compute divisive hierarchical clustering FOR PCA FEATURES
pca.diana <- diana(pca.features)
pca.diana$dc
pltree(pca.diana, cex = 0.6, hang = -1, main = "Dendrogram of diana")
rect.hclust(pca.diana, k = 2, border = 2:10)
pca.diana.clust <- cutree(pca.diana, k = 2)

pca.sil.diana <- silhouette(pca.diana.clust, dist(pca.features))
fviz_silhouette(pca.sil.diana, xlab="DIANA")

#fviz_cluster(list(data = pca.features, cluster = pca.diana.clust))
fviz_dend(pca.diana, k=2, cex=0.4)
fviz_dend(pca.agnes.avg, type="circular", cex=0.4, k=2,  main="Dendrogram of AGNES with Average Linkage - PCA")

# Visualization of cluster analysis results (3D scatterplot)
plot3d(pca.features$bcw.after.pca.x...1., pca.features$bcw.after.pca.x...2., pca.features$bcw.after.pca.x...3., col=pca.diana.clust, 
       pch=as.numeric(bcw.real.class.labels), size = 1, type='s', xlab="Dim_1", 
       ylab="Dim_2", zlab="Dim_3")
legend3d("topright", legend = c("malignant", "benign"), pch = as.numeric(bcw.real.class.labels), 
         col = pca.diana.clust, cex=1, inset=c(0.02))

#-------------------------------
## External Cluster Validation
#----------------------------------------------------------------------------------------------------------------------------
library(e1071)
library("fpc")

clust.results <- list(bcw.kmeans.labels,bcw.cluster.labels,bcw.agnes.avg.k2,bcw.diana.clust,bcw.clara.clust,
            pca.kmeans.labels,pca.cluster.labels,pca.agnes.avg.k2,pca.diana.clust,pca.clara.clust)
partition.agreement <- numeric(10)
j=1
for (cruster in clust.results) {
  matchClasses(table(cruster, bcw.real.class.labels), method="exact")
  part.agreement <- compareMatchedClasses(cruster, bcw.real.class.labels, method="exact")$diag
  print(part.agreement)
  partition.agreement[j] <- round(part.agreement, 3)
  j = j+1
}
pat.agree.val <- data.frame(Methods=c("KMeans", "PAM", "AGNES", "DIANA", "CLARA", "KMeans_PCA", "PAM_PCA", 
                                      "AGNES_PCA", "DIANA_PCA", "CLARA_PCA"),Partition_Agreement=partition.agreement)
ggplot(pat.agree.val, aes(x=Methods, y=Partition_Agreement, fill=Methods)) + geom_bar(stat="identity") + 
  geom_text(aes(label=Partition_Agreement), vjust=-0.3, size=3.5)


##PARTITION AGREEMENT: REAL LABELS AND MODEL LABELS COMPARED
clust.resultsss <- data.frame(bcw.kmeans.labels,bcw.cluster.labels,bcw.agnes.avg.k2,bcw.diana.clust,bcw.clara.clust,
                      pca.kmeans.labels,pca.cluster.labels,pca.agnes.avg.k2,pca.diana.clust,pca.clara.clust)
pat.res.matrix <- matrix(0, nrow = length(clust.results), ncol = length(clust.results))
colnames(pat.res.matrix) <- c("KMeans", "PAM", "AGNES", "DIANA", "CLARA", "KMeans_PCA", "PAM_PCA", 
                                      "AGNES_PCA", "DIANA_PCA", "CLARA_PCA")
rownames(pat.res.matrix) <- c("KMeans", "PAM", "AGNES", "DIANA", "CLARA", "KMeans_PCA", "PAM_PCA", 
                                      "AGNES_PCA", "DIANA_PCA", "CLARA_PCA")
pat.resprand.matrix <- matrix(0, nrow = length(clust.results), ncol = 1) #FOR THE RAND INDEX
acc.res.vector <- matrix(0, nrow = length(clust.results), ncol = 1) #ACCURACY
species <- as.numeric(bcw.real.class.labels)

for (i in 1:length(clust.results)){
  for (j in 1:length(clust.results)){
    if (i==j){
      part.agreement <- compareMatchedClasses(clust.resultsss[,i], bcw.real.class.labels, method="exact")$diag
      pat.res.matrix[i,j] <- round(part.agreement, 3)
      accuracy <- mean(species == clust.resultsss[,i])
      if (i==1){
        acc.res.vector[i] <- round(1-accuracy, 3)
      } else{
        acc.res.vector[i] <- round(accuracy, 3) 
      }
      if (i<6){
        clust_stats <- cluster.stats(d = dist(pca.features), species, clust.resultsss[,i])
        pat.resprand.matrix[i] <- round(clust_stats$corrected.rand, 3)
      } else {
        clust_stats <- cluster.stats(d = dist(bcw.features), species, clust.resultsss[,i])
        pat.resprand.matrix[i] <- round(clust_stats$corrected.rand, 3)
      }
    } else {
      part.agreement <- compareMatchedClasses(clust.resultsss[,i], clust.resultsss[,j], method="exact")$diag
      pat.res.matrix[i,j] <- round(part.agreement, 3)
    }
  }
}
library(pheatmap)
pheatmap(pat.res.matrix, display_numbers = T, cluster_rows = F, cluster_cols = F, number_format = "%.3f")
pat.res.matrix #RAND INDEX
acc.res.vector #ACCURACY MEASURE






