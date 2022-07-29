####S2A####
DimPlot(combined, split.by = "orig.ident", reduction = "tsne")
####S2B####
as.data.frame(table(combined@meta.data$final_clusters, combined@meta.data$kit, combined@meta.data$orig.ident))
