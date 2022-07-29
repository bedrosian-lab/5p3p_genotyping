#### S1 ####
#Load Integrated Object
combined <- readRDS("~/02_MethodsProject/12_FinalRScripts/combined.aftcelltype.rds")

##S1A########
#Total # of nuclei per sample
#Made data frame based on cell ranger info
as.data.frame(table(combined@meta.data$orig.ident, combined@meta.data$kit))
totalnumberofnucleipsdf <- data.frame(Patient_Population = rep(c("47", "48", "49"), each = 2),
                                      Count=c(6223, 5864, 8858, 8000, 5783, 5971),
                                      ChromiumKit=rep(c("3p", "5p"), 3))
head(totalnumberofnucleipsdf)

#Made plot
totalnumberofnucleipsplot <- ggplot(totalnumberofnucleipsdf, aes(x =Patient_Population, y =Count, fill =ChromiumKit)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("black", "grey")) + ggtitle("Total Number of Nuclei")
totalnumberofnucleipsplot


##S1B########
#Genes detected
#Made data frame based on cell ranger info
genesdetecteddf <- data.frame(Patient_Population = rep(c("47", "48", "49"), each = 2),
                              Count=c(32723, 31480, 31058, 29424, 31167, 30913),
                              ChromiumKit=rep(c("3p", "5p"), 3)) + NoLegend()


#Made plot 
Genedetectedplot <- ggplot(genesdetecteddf, aes(x =Patient_Population, y =Count, fill =ChromiumKit)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("black", "grey")) + ggtitle("Total Genes Detected") 

##S1C########
#Total reads 
#Made data frame based on cell ranger info
readsdetecteddf <- data.frame(Patient_Population = rep(c("47", "48", "49"), each = 2),
                              Count=c(405477358, 434407783, 426072445, 366154055, 187705461, 286175257),
                              ChromiumKit=rep(c("3p", "5p"), 3)) 


#Made plot 
readsdetectedplot <- ggplot(genesdetecteddf, aes(x =Patient_Population, y =Count, fill =ChromiumKit)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("black", "grey")) + ggtitle("Total Reads") 
readsdetectedplot

#### S1D -S1F ####

##### S1D ####
#47
p1 <- VlnPlot(p47, group.by = "orig.ident", features = "nCount_RNA", pt.size = 0) + labs(x = "Method") + NoLegend()
p1
p2 <- VlnPlot(p47, group.by = "orig.ident", features = "nFeature_RNA", pt.size = 0) + labs(x = "Method") + NoLegend()
p3 <- VlnPlot(p47, group.by = "orig.ident", features = "percent.mt", pt.size = 0) + labs(x = "Method") + NoLegend()
p4 <- VlnPlot(p47, group.by = "orig.ident", features = "percent.ribo", pt.size = 0) + labs(x = "Method")

#combined sensitivity plots
sensitivityplot47 <- plot_grid(p1, p2, p3,p4, ncol = 4)
sensitivityplot47

##### S1E ####
#48
p5 <- VlnPlot(p48, group.by = "orig.ident", features = "nCount_RNA", pt.size = 0) + labs(x = "Method") + NoLegend()
p6 <- VlnPlot(p48, group.by = "orig.ident", features = "nFeature_RNA", pt.size = 0) + labs(x = "Method") + NoLegend()
p7 <- VlnPlot(p48, group.by = "orig.ident", features = "percent.mt", pt.size = 0) + labs(x = "Method") + NoLegend()
p8 <- VlnPlot(p48, group.by = "orig.ident", features = "percent.ribo", pt.size = 0) + labs(x = "Method") 

#combined sensitivity plots 
sensitivityplot48 <- plot_grid(p5,p6,p7,p8, ncol = 4)
sensitivityplot48

##### S1F ####
#49
p9 <- VlnPlot(p49, group.by = "orig.ident", features = "nCount_RNA", pt.size = 0) + labs(x = "Method") + NoLegend()
p10 <- VlnPlot(p49, group.by = "orig.ident", features = "nFeature_RNA", pt.size = 0) + labs(x = "Method") + NoLegend()
p11 <- VlnPlot(p49, group.by = "orig.ident", features = "percent.mt", pt.size = 0) + labs(x = "Method") + NoLegend()
p12 <- VlnPlot(p49, group.by = "orig.ident", features = "percent.ribo", pt.size = 0) + labs(x = "Method") + NoLegend()

#p49 combined sensitivity plots 
sensitivityplot49 <- plot_grid(p9, p10, p11, p12, ncol = 4)
sensitivityplot49

#combine all sensitivity plots
allsensitivityplot <- plot_grid(sensitivityplot47, sensitivityplot48, sensitivityplot49, ncol = 1 )
allsensitivityplot

#Arrange all plots
S1row1 <- plot_grid(totalnumberofnucleipsplot, Genedetectedplot, readsdetectedplot, ncol = 3)
S1row1
S1 <- plot_grid(S1row1,allsensitivityplot, ncol = 1)
S1
