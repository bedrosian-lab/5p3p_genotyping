####S4A####
############### CDNA LENGTH ####
#Get cDNA length for all somatic variants to make S4A
library(tidyverse)
library(GenomicFeatures)

# make 'database' from the gff we used for exome processing
txdb <- makeTxDbFromGFF("/GRCh38_latest_genomic.gff.gz")

# get transcript lengths for all tx in our txdb
txlens <- transcriptLengths(txdb) %>%
  as.data.frame() %>%
  dplyr::select(tx_name, gene_id, tx_len) %>%
  mutate(tx_gene = paste(tx_name,gene_id, sep = "_")) # make a new column called 'tx_gene' which has unique gene and transcript ID (used for joining later on)

view(txlens)
####S4B-S4D####
##Somatic Varian Analysis with Original integrated objects
#Used integrated objects from integrated analysis for each patient
#Renamed to deidentified number p47, p48, p49

#####P47####
#Load integrated object
p47 <- readRDS(file = "~/p47.combined.aftrenameIdents.rds")

####RHEB####

#Load both genotyped tables for 3p and 5p
#3p
gt3 <- read.table("~/P47_genotype_RHEB_3p.txt",
                  sep = "\t",
                  row.names = 1,
                  header = T)
view(gt3)
#5p
gt5 <- read.table("~/P47_genotype_RHEB_5p.txt",
                  sep = "\t",
                  row.names = 1,
                  header = T)
view(gt5)
#combine gt3 and gt5

#make 2 data frames (1 each SD)  where barcodes are also the first column
#3p 
df3p <- tibble::rownames_to_column(gt3, "barcodes")
#5p
df5p <- tibble::rownames_to_column(gt5, "barcodes")

#Add genotyped calls to P47 Dataset
#Extract Barcodes
metadata <- p47@assays$RNA@counts@Dimnames[[2]]

#convert from chr to data frame
metadata <- as.data.frame(metadata)

#Rename the column to metadata - easier to reference later
colnames(metadata)[1] <- "metadata"

#Combine df3p and 5p


#merging did not work as _1 is added to integrated objects as there may be similar barcodes- to fix this use the gsub command to add or remove _1 
#https://statisticsglobe.com/insert-character-pattern-in-string-r
df3p$barcodes <- gsub("^(.{18})(.*)$", "\\1_1\\2", df3p$barcodes)
df5p$barcodes <- gsub("^(.{18})(.*)$", "\\1_2\\2", df5p$barcodes)
view(df5p)

#combine 3p and 5p
#dfp47combined <- left_join(df3p, df5p, by = "barcodes")
df5_3 <- dplyr::full_join(df3p, df5p)
#left_join allows for merging to occur between furthest left columns of 2 dataframes
merged2 <- left_join(metadata, df5_3, by = c("metadata" = "barcodes"))

rownames(merged2) <- merged2$metadata

#add to objects metadata
p47 <- AddMetaData(object = p47, 
                   metadata = merged2) 
view(p47)

#Change alt/ref to alt/alt
p47$RHEB <- recode(p47$RHEB, 'alt/ref' = "alt/alt")

p47@meta.data$orig.ident <- recode(p47@meta.data$orig.ident, 'p473p' = '3p', 'p4785p' = '5p')
p47_rheb_combined_tsne <- TSNEPlot(object = p47, split.by = "orig.ident",
                                   label = T, 
                                   cols = c("azure2", "yellow","red", "black"),
                                   pt.size = 2, reduction= "tsene",
                                   group.by = "RHEB", 
                                   label.size = 0.0, 
                                   order = c("alt/alt" ,"alt/ref", "ref/ref","No Call" )) +
  ggtitle("RHEB mutation in Patient 47 GEX") 
p47_rheb_combined_tsne
#Amount of genotyped calls per kit in p47
as.data.frame(table(p47@meta.data$RHEB, p47@meta.data$orig.ident))

#marker Plot
p47_markers <- FeaturePlot(p47,
                           features = c("MOG",
                                        # oligodendrocyte marker
                                        "P2RY12",
                                        # microglia marker
                                        "MYT1",
                                        #OPC
                                        "AQP4",
                                        # astrocyte marker
                                        "ADARB2",
                                        #MITO -catch all
                                        "LHX6",
                                        #excitatory neuron marker
                                        "SLC17A7",
                                        #CGE
                                        "SKAP1",
                                        "GRIP1"), 
                           reduction = "tsne",
                           label = F)
p47_markers

####SaveRDS with RHEB in Metadata ####
saveRDS(p47, file = "~/p47withRHEB.rds")
as.data.frame(table(p47@meta.data$RHEB, p47@meta.data$orig.ident, p47@meta.data$final_clusters))

####P48####
####SLC35A2####
p48 <- readRDS("~/02_MethodsProject/12_FinalRScripts/p48.combined.aftcelltyping.rds")

#Add SLC35A2 Genotyped calls to metadata

#Extract Barcodes and make new metadata column 
#Transpose
p48@meta.data$barcodes <- p48@assays$RNA@counts@Dimnames
p48@meta.data$barcodes <- as.character(p48@assays$RNA@counts@Dimnames)


#Load both genotyped tables for 3p and 5p
#3p
gt3 <- read.table(file = "~/P48_genotype_SLC35A2_3p.txt",
                  sep = "\t",
                  header = T)
view(gt3)
#5p
gt5 <- read.table(file = "~/P48_genotype_SLC35A2_5p.txt",
                  sep = "\t",
                  header = T)
view(gt5)
#combine gt3 and gt5

#make 2 data frames (1 each SD)  where barcodes are the first column
#3p 
df3p <- tibble::rownames_to_column(gt3, "barcodes")
#5p
df5p <- tibble::rownames_to_column(gt5, "barcodes")

#Add genotyped calls to p48 Dataset
#Extract Barcodes
metadata <- p48@assays$RNA@counts@Dimnames[[2]]

#convert from chr to data frame
metadata <- as.data.frame(metadata)

#Rename the column to metadata - easier to reference later
colnames(metadata)[1] <- "metadata"

#Combine df3p and 5p
df3p$barcodes <- gsub("^(.{18})(.*)$", "\\1_1\\2", df3p$barcodes)
df5p$barcodes <- gsub("^(.{18})(.*)$", "\\1_2\\2", df5p$barcodes)
view(df5p)

df5_3 <- dplyr::full_join(df3p, df5p)

merged2 <- left_join(metadata, df5_3, by = c("metadata" = "barcodes"))

rownames(merged2) <- merged2$metadata

#add to objects metadata
p48 <- AddMetaData(object = p48, 
                   metadata = merged2) 
view(p48)

#Change alt/ref to alt/alt
p48$SLC35A2 <- recode(p48$SLC35A2, 'alt/ref' = "alt/alt")

p48_slc35a2_combined_tsne <- TSNEPlot(object = p48, split.by = "orig.ident",
                                      label = T, 
                                      cols = c("azure2", "yellow","red", "black"),
                                      pt.size = 2, 
                                      group.by = "SLC35A2", 
                                      label.size = 0.0, 
                                      order = c("alt/alt" ,"alt/ref", "ref/ref","No Call" )) +
  ggtitle("SLC35A2 mutation in Patient 48 GEX") 
p48_slc35a2_combined_tsne
#Amount of genotyped calls per kit in p48
as.data.frame(table(p48@meta.data$SLC35A2, p48@meta.data$orig.ident, p48@active.ident))
as.data.frame(table(p48@meta.data$SLC35A2, p48@meta.data$orig.ident))

#marker Plot
p48_markers <- FeaturePlot(p48,
                           features = c("MOG",
                                        # oligodendrocyte marker
                                        "P2RY12",
                                        # microglia marker
                                        "MYT1",
                                        #OPC
                                        "AQP4",
                                        # astrocyte marker
                                        "ADARB2",
                                        #MITO -catch all
                                        "LHX6",
                                        #excitatory neuron marker
                                        "SLC17A7",
                                        #CGE
                                        "SKAP1",
                                        "GRIP1"), 
                           reduction = "tsne",
                           label = F)
p48_markers

####SaveRDS with SLC35A2 in Metadata ####
saveRDS(p48, file = "~/p48withSLC35A2.rds")

####P49####
####PTPN11####
p49 <- readRDS(file = "~/02_MethodsProject/12_FinalRScripts/p49.combined.aftrenameidents.rds")
#Add PTPN11 Genotyped calls to metadata


#Extract Barcodes and make new metadata column 
#Transpose
p49@meta.data$barcodes <- p49@assays$RNA@counts@Dimnames
p49@meta.data$barcodes <- as.character(p49@assays$RNA@counts@Dimnames)


#Load both genotyped tables for 3p and 5p
#3p
gt3 <- read.table(file = "~/P49_genotype_PTPN11_3p.txt",
                  sep = "\t",
                  header = T)
view(gt3)
#5p
gt5 <- read.table(file = "~/P49_genotype_PTPN11_5p.txt",
                  sep = "\t",
                  header = T)
view(gt5)
#combine gt3 and gt5

#make 2 data frames (1 each SD)  where barcodes are the first column
#3p 
df3p <- tibble::rownames_to_column(gt3, "barcodes")
#5p
df5p <- tibble::rownames_to_column(gt5, "barcodes")

#Add genotyped calls to p49 Dataset
#Extract Barcodes
metadata <- p49@assays$RNA@counts@Dimnames[[2]]

#convert from chr to data frame
metadata <- as.data.frame(metadata)

#Rename the column to metadata - easier to reference later
colnames(metadata)[1] <- "metadata"

#Combine df3p and 5p
df3p$barcodes <- gsub("^(.{18})(.*)$", "\\1_1\\2", df3p$barcodes)
df5p$barcodes <- gsub("^(.{18})(.*)$", "\\1_2\\2", df5p$barcodes)
view(df5p)

df5_3 <- dplyr::full_join(df3p, df5p)

merged2 <- left_join(metadata, df5_3, by = c("metadata" = "barcodes"))

rownames(merged2) <- merged2$metadata

#add to objects metadata
p49 <- AddMetaData(object = p49, 
                   metadata = merged2) 
view(p49)

#Change alt/ref to alt/alt
p49$PTPN11 <- recode(p49$PTPN11, 'alt/ref' = "alt/alt")

p49_ptpn11_combined_tsne <- TSNEPlot(object = p49, split.by = "orig.ident",
                                     label = T, 
                                     cols = c("azure2", "yellow","red", "black"),
                                     pt.size = 2, 
                                     group.by = "PTPN11", 
                                     label.size = 0.0, 
                                     order = c("alt/alt" ,"alt/ref", "ref/ref","No Call" )) +
  ggtitle("PTPN11 mutation in Patient 49 GEX") 
p49_ptpn11_combined_tsne
#Amount of genotyped calls per kit in p49
as.data.frame(table(p49@meta.data$PTPN11, p49@meta.data$orig.ident, p49@active.ident))
as.data.frame(table(p49@meta.data$PTPN11, p49@meta.data$orig.ident))

#marker Plot
p48_markers <- FeaturePlot(p48,
                           features = c("MOG",
                                        # oligodendrocyte marker
                                        "P2RY12",
                                        # microglia marker
                                        "MYT1",
                                        #OPC
                                        "AQP4",
                                        # astrocyte marker
                                        "ADARB2",
                                        #MITO -catch all
                                        "LHX6",
                                        #excitatory neuron marker
                                        "SLC17A7",
                                        #CGE
                                        "SKAP1",
                                        "GRIP1"), 
                           reduction = "tsne",
                           label = F)
p48_markers

####SaveRDS with PTPN11 in Metadata ####
saveRDS(p49, file = "~/p49withPTPN11.rds")
p49 <- readRDS(file = "~/p49withPTPN11.rds")

####ROCK2####
p49 <- readRDS(file = "~/02_MethodsProject/12_FinalRScripts/p49.combined.aftrenameidents.rds")

#Chack CT labels 
DimPlot(p49, reduction = "tsne", label = F)

#Add PTPN11 Genotyped calls to metadata


#Extract Barcodes and make new metadata column 
p49@meta.data$barcodes <- p49@assays$RNA@counts@Dimnames
p49@meta.data$barcodes <- as.character(p49@assays$RNA@counts@Dimnames)

#Load both genotyped tables for 3p and 5p
#3p
gt3 <- read.table(file = "~/P49_genotype_ROCK2_3p.txt",
                  sep = "\t",
                  header = T)
view(gt3)
#5p
gt5 <- read.table(file = "~/P49_genotype_ROCK2_5p.txt",
                  sep = "\t",
                  header = T)
view(gt5)
#combine gt3 and gt5

#make 2 data frames (1 each SD)  where barcodes are the first column
#3p 
df3p <- tibble::rownames_to_column(gt3, "barcodes")
#5p
df5p <- tibble::rownames_to_column(gt5, "barcodes")

#Add genotyped calls to p49 Dataset
#Extract Barcodes
metadata <- p49@assays$RNA@counts@Dimnames[[2]]

#convert from chr to data frame
metadata <- as.data.frame(metadata)

#Rename the column to metadata - easier to reference later
colnames(metadata)[1] <- "metadata"

#Combine df3p and 5p

df3p$barcodes <- gsub("^(.{18})(.*)$", "\\1_1\\2", df3p$barcodes)
df5p$barcodes <- gsub("^(.{18})(.*)$", "\\1_2\\2", df5p$barcodes)
view(df5p)

#combine 3p and 5p
df5_3 <- dplyr::full_join(df3p, df5p)

merged2 <- left_join(metadata, df5_3, by = c("metadata" = "barcodes"))

rownames(merged2) <- merged2$metadata

#add to objects metadata
p49 <- AddMetaData(object = p49, 
                   metadata = merged2) 
view(p49)

#Change alt/ref to alt/alt
p49$ROCK2 <- recode(p49$ROCK2, 'alt/ref' = "alt/alt")

p49_rock2_combined_tsne <- TSNEPlot(object = p49, split.by = "orig.ident",
                                    label = T, 
                                    cols = c("azure2", "yellow","red", "black"),
                                    pt.size = 2, 
                                    group.by = "ROCK2", 
                                    label.size = 0.0, 
                                    order = c("alt/alt" ,"alt/ref", "ref/ref","No Call" )) +
  ggtitle("ROCK2 mutation in Patient 49 GEX") 
p49_rock2_combined_tsne
#Amount of genotyped calls per kit in P48
as.data.frame(table(p49@meta.data$ROCK2, p49@meta.data$orig.ident, p49@active.ident))
as.data.frame(table(p49@meta.data$ROCK2, p49@meta.data$orig.ident))

#marker Plot
p49_markers <- FeaturePlot(p48,
                           features = c("MOG",
                                        # oligodendrocyte marker
                                        "P2RY12",
                                        # microglia marker
                                        "MYT1",
                                        #OPC
                                        "AQP4",
                                        # astrocyte marker
                                        "ADARB2",
                                        #MITO -catch all
                                        "LHX6",
                                        #excitatory neuron marker
                                        "SLC17A7",
                                        #CGE
                                        "SKAP1",
                                        "GRIP1"), 
                           reduction = "tsne",
                           label = F)
p48_markers

####SaveRDS with ROCK2 in Metadata ####
saveRDS(p49, file = "~/p49withROCK2.rds")



