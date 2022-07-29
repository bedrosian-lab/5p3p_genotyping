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

#Load Data
#JLE493P
P483P.data <- Read10X(data.dir = "~/02_MethodsProject/01_DataFiles/JLE49_3/filtered_feature_bc_matrix/")
#JLE495P
P485P.data <- Read10X(data.dir = "~/02_MethodsProject/01_DataFiles/JLE49_5/filtered_feature_bc_matrix/")



#P48 SEURAT OBJECTS

#Create Surat Objects
#P483P
P483P <- CreateSeuratObject(counts = P483P.data, project = "p483p", min.cells = 3, min.features = 200 )
P483P
#JLE495P
P485P <- CreateSeuratObject(counts = P485P.data, project = "p485p", min.cells = 3, min.features = 200)
P485P

##Annotations
#Percent mt
P483P <- PercentageFeatureSet(P483P, pattern = "^MT-", col.name = "percent.mt")
P485P <- PercentageFeatureSet(P485P, pattern = "^MT-", col.name = "percent.mt")

#Percent Ribo
P483P <- PercentageFeatureSet(P483P, pattern = "^RP[SL]", col.name = "percent.ribo")
P485P <- PercentageFeatureSet(P485P, pattern = "^RP[SL]", col.name = "percent.ribo")

#Sample list
P483P[["sample"]] <- "P483P"
P485P[["sample"]] <- "P485P"
samples.list <- list(P483P, P485P)



#RDS files - raw
#P483P
saveRDS(P483P, file = "~/02_MethodsProject/03_AllJLE49Analysis/03_IntegratedJLE49/JLE493P.rds")
#P485P
saveRDS(P485P, file = "~/02_MethodsProject/03_AllJLE49Analysis/03_IntegratedJLE49/JLE495P.rds")
#Sample.list
saveRDS(samples.list, file = "~/02_MethodsProject/03_AllJLE49Analysis/03_IntegratedJLE49/rawsample.list.rds")

#DOUBLET FINDER
##DoubletFinder info - https://github.com/chris-mcginnis-ucsf/DoubletFinder

# compute cut-off of max Q3+1.5*IQR (boxplot max)
mt <- lapply(X = samples.list, FUN = function(x) {
  x <- x@meta.data$percent.mt
})
mt_merged <- as.data.frame(unlist(mt))
boxplot.stats(mt_merged$`unlist(mt)`)

samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- subset(x, subset = percent.mt < 5.0)
})
saveRDS(samples.list, file = "~/02_MethodsProject/02_AllJLE48Analysis/03_Integrated_JLE48_Object/JLE49samples.list.rds")

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
doubletrate.vec <- c(.067, .061)

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
saveRDS(samples.list, file = "~/02_MethodsProject/03_AllJLE49Analysis/03_IntegratedJLE49/49samples.list.withdoublets.rds")

# remove doublets from objects - MOST IMPORTANT PART
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- subset(x, subset = DF == "Singlet")
})


#RDS file- rmv Doublets
saveRDS(samples.list, file = "~/02_MethodsProject/03_AllJLE49Analysis/03_IntegratedJLE49/49samples.list.rmvdoublets.rds")

######################################################## INtegration Following DOublet Finder ##################################################

#INTEGRATION
features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list, anchor.features = features)

## Find anchors and integrate
samples.anchors <- FindIntegrationAnchors(object.list = samples.list, normalization.method = "SCT", anchor.features = features)
p48 <- IntegrateData(anchorset = samples.anchors, normalization.method = "SCT")
p48

#Save RDS file - Integrated
saveRDS(p48, file = "~/02_MethodsProject/03_AllJLE49Analysis/03_IntegratedJLE49/p48.combined.rds")
p48 <- readRDS("~/02_MethodsProject/03_AllJLE49Analysis/03_IntegratedJLE49/p48.combined.rds")
######################################### STANDARD WORKFLOW AFTER DOUBLET FINDER ################################################################################
p48 <- RunPCA(p48, verbose = T)
p48 <- RunUMAP(p48, dims = 1:30, verbose = T)
p48 <- RunTSNE(p48, dims = 1:30, verbose = T, perplexity = 30)

######################################## CELL CLUSTERING AFTER DOUBLET FINDER ##################################
p48 <- FindNeighbors(p48, reduction = "pca", dims = 1:20)
p48 <- FindClusters(p48, resolution = 0.5)

#Find Markers
p48.markers <- FindAllMarkers(p48, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
p48.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


#RENAME IDENTS
Idents(object = p48) <- "seurat_clusters"
#Check active Ident in console
p48@active.ident

#Rename Idents
p48 <- RenameIdents(p48,'0' = 'Oligodendrocytes', '1' = 'Oligodendrocytes', '2' = 'Oligodendrocytes', '3' = 'Oligodendrocytes', '4' = 'Oligodendrocytes', '5' = 'Oligodendrocytes', '6' = 'Microglia', '7' = 'OPCs', '8' = 'Oligodendrocytes', '9' = 'Microglia', '10' = 'Oligodendrocytes', '11' = 'Astrocytes', '12' = 'Oligodendrocytes', '13' = 'CGE-derived Interneurons', '14' = 'Excitatory Neurons', '15' = 'OPCs')
#Confirm cluster reidentifiction
DimPlot(p48, split.by = "sample", reduction = "umap", label = T) 

#In this case split MGE and Lymphocytes
plot <- DimPlot(p48, reduction = "umap", label = FALSE)
HoverLocator(plot = plot, information = FetchData(p48, vars = c("ident", "PC_1", "nFeature_RNA")))
#Allows you to select a population of cells
select.cells <- CellSelector(plot = plot)  

Idents(p48, cells = select.cells) <- "MGE-derived Interneurons"


plot <- DimPlot(p48, reduction = "umap", label = FALSE)
HoverLocator(plot = plot, information = FetchData(p48, vars = c("ident", "PC_1", "nFeature_RNA")))
#Allows you to select a population of cells
select.cells <- CellSelector(plot = plot)  

Idents(p48, cells = select.cells) <- "Lymphocytes"

#Add Labels to metadata
p48$final_clusters <- Idents(p48)

#RDS- After CTyping
saveRDS(p48, file = "~/02_MethodsProject/12_FinalRScripts/p48.combined.aftcelltyping.rds")




