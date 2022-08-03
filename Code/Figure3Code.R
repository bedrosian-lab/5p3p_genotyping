#Load Integrated Object
combined <- readRDS("~/02_MethodsProject/12_FinalRScripts/combined.aftcelltype.rds")
#Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(viridis)

####3A####
DimPlot(combined, reduction = "tsne")

####3B####
DimPlot(combined, reduction = "tsne", group.by = "kit")

####3C####
#Generate feature plots for each marker gene
FeaturePlot(combined, reduction = "tsne", features = "MOG", cols = c("grey", "#F8766D"))
FeaturePlot(combined, reduction = "tsne", features = "P2RY12", cols = c("grey", "goldenrod"))
FeaturePlot(combined, reduction = "tsne", features = "MYT1", cols = c("grey", "#7CAE00"))
FeaturePlot(combined, reduction = "tsne", features = "AQP4", cols = c("grey", "#00b159"))
FeaturePlot(combined, reduction = "tsne", features = "ADARB2", cols = c("grey", "#00A5FF"))
FeaturePlot(combined, reduction = "tsne", features = "LHX6", cols = c("grey", "#FC61D5"))
FeaturePlot(combined, reduction = "tsne", features = "SKAP1", cols = c("grey", "#DC71FA"))
FeaturePlot(combined, reduction = "tsne", features = "SLC17A7", cols = c("grey", "#00aBf7"))

####3D####
av.exp <- AverageExpression(combined, assays = "RNA", slot = "data", group.by = c("kit", "final_clusters"))
av.exp.df <- as.data.frame(av.exp$RNA)
cor.exp <- as.data.frame(cor(av.exp.df))
pheatmap(cor.exp, color=colorRampPalette(c("#D44292", "#952AE0", "#4B2991"))(50), angle_col = 90, fontsize = 14, cluster_rows = F, cluster_cols = F) 

####3E####
#Get Cell Type Proportions for plot made in Graphpad Prism
as.data.frame(table(combined@meta.data$final_clusters, combined@meta.data$kit))

####3F####

############## MARKER SCORE #################
# Get the top 10 marker genes for each cluster
top10 <- read.csv("~/02_MethodsProject/10_MarkerScore/top10.csv")

# Get CPMs per https://satijalab.org/seurat/reference/normalizedata, They will be stored in RNA assay, counts slot.
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined, normalization.method = "RC", scale.factor = 1e6)

# Make a list of the marker genes to use for calculating proportions
genes <- unique(top10$gene)

result <- c()
clusters <- levels(combined@meta.data$final_clusters)

#Subset 3P and 5P - So marker score comparison can be between 3p and 5p 

#subset by 3'
all3 <- subset(combined, subset = kit == "3p")
head(all3)
#subset by 5'
all5 <- subset(combined, subset = kit == "5p")
head(all5)

#Marker Score Loop 3'
for (i in 1:length(clusters)){
  subset <- subset(all3, final_clusters == clusters[i])
  DefaultAssay(subset) <- "RNA"
  prop_cluster <- c()
  for (i in 1:length(genes)){
    gene <- FetchData(subset, vars = genes[i], slot = "counts")
    prop <- (sum(gene[[1]] >1)) / (nrow(gene))
    prop_cluster <- append(prop_cluster, prop)
  }
  result <- cbind(result, prop_cluster)
}
result <- as.data.frame(result)
rownames(result) <- genes
colnames(result) <- clusters
# result now has the proportions of cells from the 3' kit expressing each marker gene >1 CPM for each cluster
write.csv(result, file = "~/02_MethodsProject/10_MarkerScore/3pmarkerscore.csv")

#Marker score loop 5' 
for (i in 1:length(clusters)){
  subset <- subset(all5, final_clusters == clusters[i])
  DefaultAssay(subset) <- "RNA"
  prop_cluster <- c()
  for (i in 1:length(genes)){
    gene <- FetchData(subset, vars = genes[i], slot = "counts")
    prop <- (sum(gene[[1]] >1)) / (nrow(gene))
    prop_cluster <- append(prop_cluster, prop)
  }
  result <- cbind(result, prop_cluster)
}
result <- as.data.frame(result)
rownames(result) <- genes
colnames(result) <- clusters
# result now has the proportions of cells from the 5' kit expressing each marker gene >1 CPM for each cluster
write.csv(result, file = "~/02_MethodsProject/10_MarkerScore/5pmarkerscore.csv")

#After Determination of cells expressing above background (>1 CPM) marker score calculations were performed in excel
#After determination of Marker Score for each files read in to find correlation between capture methods

#Read in csv of calculations 
markerscore3p <- read.csv("~/02_MethodsProject/10_MarkerScore/markerscore3p.csv")
head(markerscore3p)

markerscore5p <- read.csv("~/02_MethodsProject/10_MarkerScore/markerscore5p.csv")
head(markerscore5p)

#Make markerscore correlation plot
#Create dataframe with 3p marker score, 5p marker score, and labels
markerscoredf <- data.frame(x=markerscore3p$Sum.2.Sum, y=markerscore5p$Sum.2.Sum, z=markerscore3p$Gene)
#
ggplot(markerscoredf, aes(x,y)) +
  geom_point() +
  geom_text_repel(aes(label = z)) +
  title(main = "3p vs. 5' Marker Score", xlab = "3p GEX Kit Marker Score", ylab = "5' GEX kit Marker Score")

ggplot(markerscoredf, aes(x, y)) +
  geom_point() +
  geom_text_repel(aes(label = z)) +
  ggtitle("3p vs. 5p Marker Score") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("3p GEX Kit Marker Score") + 
  ylab("5p GEX kit Marker Score") + theme(axis.title=element_text(size=16), plot.title = element_text(size=18))
# get Rsquared value for correlation, can add to plot manually
#Join Data Tables
join <- dplyr::full_join(markerscore3p, markerscore5p)
summary(lm(markerscore3p$Sum.2.Sum ~ markerscore5p$Sum.2.Sum, data = join))


