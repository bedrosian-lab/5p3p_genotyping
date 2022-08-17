#Dependencies
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


#P49 5&3 P Objects

##Load Data
P493P.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
P495P.data <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")

##Create Objects
P493P <- CreateSeuratObject(P493P.data, project = "p493p", min.cells = 3, min.features = 200)
P495P <- CreateSeuratObject(P495P.data, project = "p495p", min.cells = 3, min.features = 200)


##Annotations
#perceent mt
P493P <- PercentageFeatureSet(P493P, pattern = "^MT-",  col.name = "percent.mt")
P495P <- PercentageFeatureSet(P495P, pattern = "^MT-",  col.name = "percent.mt")
#percent.ribo
P493P <- PercentageFeatureSet(P493P, pattern = "^RP[SL]", col.name = "percent.ribo")
P495P <- PercentageFeatureSet(P495P, pattern = "^RP[SL]", col.name = "percent.ribo")

#sample list
P493P[["sample"]] <- "P493P"
P495P[["sample"]] <- "P495P"
samples.list <- list(P493P, P495P)


#DOUBLET FINDER
##DoubletFinder info - https://github.com/chris-mcginnis-ucsf/DoubletFinder

samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- subset(x, subset = percent.mt < 5)
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
doubletrate.vec <- c(.048, .048)

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

#remove doublets from objects
samples.list <- lapply(X = samples.list, FUN = function(x) {
x <- subset(x, subset = DF == "Singlet")
})

#Plot after Doublet removal
for (i in 1:length(samples.list)) {
  plot2 <- DimPlot(samples.list[[i]], reduction = "umap", group.by = "DF", pt.size = .7)
}
plot2 

#PLot with and without Doublet
plot + plot2

######################################################## INtegration Following DOublet Finder ##################################################

features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list, anchor.features = features)

## Find anchors and integrate
samples.anchors <- FindIntegrationAnchors(object.list = samples.list, normalization.method = "SCT", anchor.features = features)
p49 <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT")
p49


######################################### STANDARD WORKFLOW AFTER DOUBLET FINDER ################################################################################
p49 <- RunPCA(p49, verbose = T)
p49 <- RunUMAP(p49, dims = 1:30, verbose = T)
p49 <- RunTSNE(p49, dims = 1:30, verbose = T, perplexity = 30)



######################################## CELL CLUSTERING AFTER DOUBLET FINDER ##################################
p49 <- FindNeighbors(p49, reduction = "pca", dims = 1:20)
p49 <- FindClusters(p49, resolution = 0.5)

#Find Markers
p49.markers <- FindAllMarkers(p49, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
p49.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#RENAME IDENTS
#Check active Ident in console
p49@active.ident

#Rename Idents
p49 <- RenameIdents(p49, '0' = "Microglia", '1' = "OPCs", '2' = "Astrocytes", '3' = "Oligodendrocytes", '4' = "Excitatory Neurons", '5' = "Oligodendrocytes", '6' = "Microglia", '7' = "Excitatory Neurons", '8' = "Mitochondrial", '9' = "MGE-derived Interneurons", '10' = "Excitatory Neurons", '11' = "OPCs", '12' = "MGE-derived Interneurons", '13' = "CGE-derived Interneurons", '14' = "Excitatory Neurons", '15' = "Mitochondrial", '16' = "OPCs")

plot <- DimPlot(p49, reduction = "umap", label = T)
HoverLocator(plot = plot, information = FetchData(p49, vars = c("ident", "PC_1", "nFeature_RNA")))
#Allows you to select a population of cells
select.cells <- CellSelector(plot = plot)  

Idents(p49, cells = select.cells) <- "Lymphocytes"

#Add Labels to metadata
p49$final_clusters <- Idents(p49)

##Dimplot new Idents
DimPlot(p49, split.by = "orig.ident", repel = TRUE, label.size = 3, reduction = "tsne")

#RDS- After RenameIdents
saveRDS(p49, file = "~/02_MethodsProject/12_FinalRScripts/p49.combined.aftrenameidents.rds")

