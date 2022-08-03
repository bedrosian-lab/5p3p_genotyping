library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(viridis)

#Load Data
combined <- readRDS("~/02_MethodsProject/05_AllSampleAnalysis/combined.aftcelltype.rds")
####2A####

#Cell Ranger Outputs

####2B####
VlnPlot(combined,
        features = "nFeature_RNA",
        group.by = "kit",
        pt.size = 0) + theme(axis.text.x = element_text(angle = 0) + NoLegend()) +  theme(plot.title = element_blank()) + labs(x = "Method", y = "UMIs detected per nucleus")

####2C####
VlnPlot(combined,
        features = "nCount_RNA",
        group.by = "kit",
        pt.size = 0) + theme(axis.text.x = element_text(angle = 0)) + NoLegend() +  theme(plot.title = element_blank()) + labs(x = "Method", y = "Genes detected per nucleus")


####2D####
VlnPlot(combined,
          features = "percent.mt",
          group.by = "kit",
          pt.size = 0) + theme(axis.text.x = element_text(angle = 0)) + NoLegend() + theme(plot.title = element_blank()) + labs(x = "Method", y = "Percentage of mitochondrial reads per nuc.")
            
####2E####
VlnPlot(combined,
        features = "percent.ribo",
        group.by = "kit",
        pt.size = 0) + theme(axis.text.x = element_text(angle = 0)) + NoLegend() + theme(plot.title = element_blank()) + labs(x = "Method", y = "Percentage of ribosomal reads per nuc.")


####2F####
i3p <- subset(combined, subset = kit=="3p")
DefaultAssay(i3p) <- "RNA"
i3p <- NormalizeData(i3p)
Avg <- AverageExpression(i3p, assays = "RNA", slot = "data", group.by = "kit")
rownames(Avg$RNA) -> Gene
df <- as.data.frame(unlist(Avg$RNA))
df$Gene <- Gene

i5p <- subset(combined, subset = kit=="5p")
DefaultAssay(i5p) <- "RNA"
i5p <- NormalizeData(i5p)
Avg <- AverageExpression(i5p, assays = "RNA", slot = "data", group.by = "kit")
rownames(Avg$RNA) -> Gene
df_5p <- as.data.frame(unlist(Avg$RNA))
df_5p$Gene <- Gene

join <- inner_join(df, df_5p, by = "Gene")
ggplot(join, aes(x= all.x, y=all.y, color="black")) + geom_point(size=1) + xlim(0,60) + ylim(0,60) + xlab("Mean log normalized expression in 3p assay") + ylab("Mean log normalized expression in 5p assay") + theme_bw() + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + NoLegend()

# get Rsquared value for correlation, can add to plot manually, R2=0.87
summary(lm(all.y ~ all.x, data = join))

####2G####
#Patient Specific Correlations
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
av.exp <- AverageExpression(combined, assays = "RNA", slot = "data", group.by = "orig.ident")
av.exp.df <- as.data.frame(av.exp$RNA)

#47
ggplot(av.exp.df, aes(x=p47_3p, y=p47_5p, color="black")) + geom_point(size=1) + xlim(0,60) + ylim(0,60) + theme_bw() + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                                                                                                                                                              axis.title.y = element_blank()) + labs(title="Patient 47") + theme(plot.title = element_text(hjust = 0.5))
#Get Rsquared Value for patient 47 
summary(lm(av.exp.df$`47_3p` ~ av.exp.df$'47_5p', data = av.exp.df))
#r^2 = .91


#48
ggplot(av.exp.df, aes(x=p48_3p, y=p48_5p, color="black")) + geom_point(size=1) + xlim(0,60) + ylim(0,60) + theme_bw() + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                                                                                                                                                              axis.title.y = element_blank()) + labs(title="Patient 48") + theme(plot.title = element_text(hjust = 0.5))
#Get Rsquared Value for patient 47 
#Get Rsquared Value for patient 47 
summary(lm(av.exp.df$p48_3p ~ av.exp.df$p48_5p, data = av.exp.df))
#r^2 = .97

#49
ggplot(av.exp.df, aes(x=p49_3p, y=p49_5p, color="black")) + geom_point(size=1) + xlim(0,60) + ylim(0,60) + theme_bw() + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                                                                                                                                                              axis.title.y = element_blank()) + labs(title="Patient 49") + theme(plot.title = element_text(hjust = 0.5))
#Get Rsquared Value for patient 47 
#Get Rsquared Value for patient 47 
summary(lm(av.exp.df$p49_3p ~ av.exp.df$p49_5p, data = av.exp.df))
#r^2 = .95

####2H####
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
av.exp <- AverageExpression(combined, assays = "RNA", slot = "data", group.by = "orig.ident")
av.exp.df <- as.data.frame(av.exp$RNA)
cor.exp <- as.data.frame(cor(av.exp.df))
pheatmap(cor.exp, color=viridis(50), angle_col = 0, fontsize = 14)

