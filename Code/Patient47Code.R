library(Seurat)
library(dplyr)
library(patchwork)
library(glmGamPoi)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(Matrix)
library(parallel)
library(fields)
library(cowplot)
library(scCATCH)
library(KernSmooth)
library(ROCR)
library(sctransform)

######## p47 3p vs. 5p Integrated Object #########

#Load Datasets#
#Load p47_3P Data
p47_3P.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
#Create Object
P47_3P <- CreateSeuratObject(counts = p47_3P.data, project = "p473p", min.cells = 3, min.features = 200)
P47_3P

#Load p47_5P Data
p47_5P.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix")
#Create Object
P47_5P <- CreateSeuratObject(counts = p47_5P.data, project = "p475p", min.cells = 3, min.features = 200)
P47_5P

#Add annotations

#Percent Mitochondrial
P47_3P <- PercentageFeatureSet(P47_3P, pattern = "^MT-", col.name = "percent.mt")
P47_5P <- PercentageFeatureSet(P47_5P, pattern = "^MT-", col.name = "percent.mt")
#Percent Ribosomal
P47_3P <- PercentageFeatureSet(P47_3P, pattern = "^RP[SL]", col.name = "percent.ribo")
P47_5P <- PercentageFeatureSet(P47_5P, pattern = "^RP[SL]", col.name = "percent.ribo")

#DOUBLET FINDER
#DoubletFinder info - https://github.com/chris-mcginnis-ucsf/DoubletFinder

#Make Sample list
P47_3P[["sample"]] <- "3p"
P47_5P[["sample"]] <- "5p"
samples.list <- list(P47_3P, P47_5P)

samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- subset(x, subset = percent.mt < 5.0)
})
#Standard QC workflow
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt")
  x <- RunPCA(x, verbose = T)
  x <- RunUMAP(x, dims = 1:30)
  x <- RunTSNE(x, dims = 1:30)
})

# run sweep stats
sweep.stats.list <- list()
pk.vec <- vector()
for (i in 1:length(samples.list)) {
  samples_temp <- samples.list[[i]]
  sweep.res.list <- paramSweep_v3(samples_temp, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk.vec[[i]] <- as.numeric(as.vector(top_n(bcmvn,1,BCmetric)[["pK"]]))
  sweep.stats.list[[i]] <- sweep.stats
}

# manually put in your expected doublet rate for each sample based on 10x genomics documentation
doubletrate.vec <- c(.056, .046)

# run doublet finder for each sample
for (i in 1:length(samples.list)) {
  samples_temp <- samples.list[[i]]
  nExp_poi <- doubletrate.vec[i]*nrow(samples_temp@meta.data)
  samples_temp <- doubletFinder_v3(samples_temp, PCs = 1:30, pN = 0.25, pK = pk.vec[i], nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  samples.list[[i]] <- samples_temp
}

# create a metadata category named "DF" in each Seurat object that is the same as the DF.classifications_... 
# that will combine when the sets are integrated
#Makes if easier to call later
for (i in 1:length(samples.list)) {
  samples.list[[i]]@meta.data[["DF"]] <- samples.list[[i]]@meta.data[[10]]
}

# manual step - take a look at the doublets on a plot
for (i in 1:length(samples.list)) {
  plot <- DimPlot(samples.list[[i]], reduction = "umap", group.by = "DF", pt.size = .7)
}
plot
plot + NoLegend()

# remove doublets from objects - MOST IMPORTANT PART
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- subset(x, subset = DF == "Singlet")
})

######################################################## Integration Following Doublet Finder ##################################################
## Finish prepping for integration the SCT way
features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list, anchor.features = features)

## Find anchors and integrate
samples.anchors <- FindIntegrationAnchors(object.list = samples.list, normalization.method = "SCT", anchor.features = features)
p47 <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT")


######################################### STANDARD WORKFLOW AFTER DOUBLET FINDER ################################################################################
p47 <- RunPCA(p47, verbose = T)
p47 <- RunUMAP(p47, dims = 1:30, verbose = T)
p47 <- RunTSNE(p47, dims = 1:30, verbose = T, perplexity = 30)

######################################## CELL CLUSTERING AFTER DOUBLET FINDER ##################################
p47 <- FindNeighbors(p47, reduction = "pca", dims = 1:20)
p47 <- FindClusters(p47, resolution = 0.5)

#Find Markers
p47.markers <- FindAllMarkers(p47, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
p47.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# CELL TYPING
#Duplicate seurat clusters
p47@meta.data$seurat_clusters2 <- p47@meta.data$seurat_clusters

#Change active Ident to seurat_clusters2
Idents(object= p47) <- "seurat_clusters2"

#Confirm active Ident Change
p47@active.ident

# Identify the clusters as cell types
p47 <- RenameIdents(p47,'0' = 'Microglia', '1' = 'OPC', '2' = 'Oligodendrocytes', '3' = 'Excitatory Neurons', '4' = 'Excitatory Neurons', '5' = 'Excitatory Neurons', '6' = 'CGE-derived Interneurons', '7' = 'Lymphocytes', '8' = 'MGE-derive Interneurons', '9' = 'Mitochondrial', '10' = 'OPC', '11' = 'Astrocytes', '12' = 'OPC', '13' = 'Excitatory Neurons', '14' = 'Excitatory Neurons', '15' = 'CGE-derived Imternuerons', '16' = 'MGE-derived Interneurons', '17' = 'Lymphocytes')

#Confirm cluster reidentifiction
DimPlot(p47, split.by = "sample") 

#Add Labels to metadata
p47$final_clusters <- Idents(p47)
#RDS after RenameIdents

saveRDS(p47, file = "~/02_MethodsProject/12_FinalRScripts/p47.combined.aftrenameIdents.rds")

