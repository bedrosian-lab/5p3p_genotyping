#Load Combined object
combined <- readRDS("~/02_MethodsProject/12_FinalRScripts/combined_all_nodoubletsaftint.rds")

####4A####
######### GERMLINE VARIANT ########
# # ## # # # # # # # ## # #. # # # # #. ##. # # # # # # # # 
# Plotting percent of cells genotyped by position for all germline variants
# # ## # # # # ## # # # # ## # #. ## # # # # # # # # # # # # 
#For 3’ you have to subtract the transcript length minus the nucleotide position of the variant:
p47_germline_3p_geno <- read.table("/igm/projects/Singlecell_analysis_projects/JLE_snRNA_objects/Tracy_Katie_JLE_snRNA/JLE48_germline_3p_geno.txt", header=T, sep="\t")
p48_germline_3p_geno <- read.table("/igm/projects/Singlecell_analysis_projects/JLE_snRNA_objects/Tracy_Katie_JLE_snRNA/JLE49_germline_3p_geno.txt", header=T, sep="\t")
p49_germline_3p_geno <- read.table("/igm/projects/Singlecell_analysis_projects/JLE_snRNA_objects/Tracy_Katie_JLE_snRNA/JLE50_germline_3p_geno.txt", header=T, sep="\t")

join <- rbind(p47_germline_3p_geno, p48_germline_3p_geno)
join <- rbind(join, p49_germline_3p_geno)

join$dist <- join$tx_len - join$cDNA_HGVS
ggplot(join, aes(x=dist, y=perc_geno)) + geom_point(size=1) + xlab("Distance of variant from 3' transcript end") + ylab("Percent of cells genotyped") + theme_bw()
# might show a zoomed up version to make the point that we can genotype outside the '90' bases you're supposed to sequence. maybe because of non specific poly A capture sites
#ggplot(join, aes(x=dist, y=perc_geno)) + geom_point(size=1) + xlab("Distance of variant from 3' transcript end") + ylab("Percent of cells genotyped") + xlim(0,6000) + theme_bw()
# we can pull out an example in IGV to show the coverage far from the transcript end


#For 5’ you can just use the nucleotide position
#I'm going to join all three datasets because I think it's more important to show a larger variant set than patient to patient differences (none are expected)
p47_germline_5p_geno <- read.table("/igm/projects/Singlecell_analysis_projects/JLE_snRNA_objects/Tracy_Katie_JLE_snRNA/JLE48_germline_5p_geno.txt", header=T, sep="\t")
p48_germline_5p_geno <- read.table("/igm/projects/Singlecell_analysis_projects/JLE_snRNA_objects/Tracy_Katie_JLE_snRNA/JLE49_germline_5p_geno.txt", header=T, sep="\t")
p49_germline_5p_geno <- read.table("/igm/projects/Singlecell_analysis_projects/JLE_snRNA_objects/Tracy_Katie_JLE_snRNA/JLE50_germline_5p_geno.txt", header=T, sep="\t")
join <- rbind(p47_germline_5p_geno, p48_germline_5p_geno)
join <- rbind(join, p49_germline_5p_geno)
ggplot(join, aes(x=cDNA_HGVS, y=perc_geno)) + geom_point(size=1) + xlab("Distance of variant from 5' transcript end") + ylab("Percent of cells genotyped") + theme_bw()
# might show a zoomed up version to make the point that we can genotype outside the '90' bases you're supposed to sequence. maybe because of non specific poly A capture sites
ggplot(join, aes(x=cDNA_HGVS, y=perc_geno)) + geom_point(size=1) + xlab("Distance of variant from 5' transcript end") + ylab("Percent of cells genotyped") + xlim(0,1000) + theme_bw()
# we can pull out an example in IGV to show the coverage far from the transcript end

#Writing data for Prism plots
p47 <- inner_join(p47_germline_3p_geno, p47_germline_5p_geno, by = "variant")
p48 <- inner_join(p48_germline_3p_geno, p48_germline_5p_geno, by = "variant")
p49 <- inner_join(p49_germline_3p_geno, p49_germline_5p_geno, by = "variant")

join <- rbind(p47, p48)
join <- rbind(join, p49)

join$dist <- join$tx_len.x - join$cDNA_HGVS.x
write.csv(join, "/igm/home/tab013/5p3p/germline.csv")

####4B####
#Plot distance by expression, with color as percentage
av.exp <- AverageExpression(combined, assays = "RNA", slot = "data", group.by = "orig.ident")
av.exp.df <- as.data.frame(av.exp$RNA)
av.exp.df$Gene <- rownames(av.exp.df)
p47_germline_3p_geno$dist <- p47_germline_3p_geno$tx_len - p47_germline_3p_geno$cDNA_HGVS
p47_3p <- inner_join(p47_germline_3p_geno, av.exp.df, by = "Gene")
p47_3p <- p47_3p[,c(1,2,11,13,14)]
p47_5p <- inner_join(p47_germline_5p_geno, av.exp.df, by = "Gene")
p47_5p <- p47_5p[,c(1,2,4,11,14)]
p48_germline_3p_geno$dist <- p48_germline_3p_geno$tx_len - p48_germline_3p_geno$cDNA_HGVS
p48_3p <- inner_join(p48_germline_3p_geno, av.exp.df, by = "Gene")
p48_3p <- p48_3p[,c(1,2,11,13,15)]
p48_5p <- inner_join(p48_germline_5p_geno, av.exp.df, by = "Gene")
p48_5p <- p48_5p[,c(1,2,4,11,16)]
p49_germline_3p_geno$dist <- p49_germline_3p_geno$tx_len - p49_germline_3p_geno$cDNA_HGVS
p49_3p <- inner_join(p49_germline_3p_geno, av.exp.df, by = "Gene")
p49_3p <- p49_3p[,c(1,2,11,13,17)]
p49_5p <- inner_join(p49_germline_5p_geno, av.exp.df, by = "Gene")
p49_5p <- p49_5p[,c(1,2,4,11,18)]

colnames(p47_3p) <- c("Variant", "Gene", "PercentGeno", "Distance", "Exp")
colnames(p47_5p) <- c("Variant", "Gene", "Distance", "PercentGeno", "Exp")
p47 <- rbind(p47_3p, p47_5p)
colnames(p48_3p) <- c("Variant", "Gene", "PercentGeno", "Distance", "Exp")
colnames(p48_5p) <- c("Variant", "Gene", "Distance", "PercentGeno", "Exp")
p48 <- rbind(p48_3p, p48_5p)
colnames(p49_3p) <- c("Variant", "Gene", "PercentGeno", "Distance", "Exp")
colnames(p49_5p) <- c("Variant", "Gene", "Distance", "PercentGeno", "Exp")
p49 <- rbind(p49_3p, p49_5p)

join <- rbind(p47, p48)
join <- rbind(join, p49)
join <- join[order(as.numeric(factor(join$PercentGeno))),]
ggplot(join, aes(x=Distance, y=Exp, color=PercentGeno)) + geom_point() + scale_colour_gradientn(colours = c("grey", "firebrick1", "firebrick2", "firebrick3", "firebrick")) + theme_classic() + xlim(0,50000) + xlab("Distance from either transcript end") + ylab("Mean log normalized expression") + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + labs(color="Percentage of Nuclei Genotyped") + theme(legend.position = "top")


####4C####
############### CDNA LENGTH ####
#Pulled out cDNA length of RHEB Variant to make plot 
library(tidyverse)
library(GenomicFeatures)

# make 'database' from the gff we used for exome processing
txdb <- makeTxDbFromGFF("/igm/projects/Singlecell_analysis_projects/JLE_snRNA_objects/Tracy_Katie_JLE_snRNA/GRCh38_latest_genomic.gff.gz")

# get transcript lengths for all tx in our txdb
txlens <- transcriptLengths(txdb) %>%
  as.data.frame() %>%
  dplyr::select(tx_name, gene_id, tx_len) %>%
  mutate(tx_gene = paste(tx_name,gene_id, sep = "_")) # make a new column called 'tx_gene' which has unique gene and transcript ID (used for joining later on)

view(txlens)

####4D####
##Somatic Varian Analysis with Original integrated objects
#Used integrated objects from integrated analysis for each patient
#Renamed to deidentified number p47, p48, p49

#####P47####
#Load integrated object
p47 <- readRDS(file = "~/02_MethodsProject/12_FinalRScripts/p47.combined.aftrenameIdents.rds")
####RHEB####

#Load both genotyped tables for 3p and 5p
#3p
gt3 <- read.table("~/02_MethodsProject/01_DataFiles/P47_3p/P47_genotype_RHEB_3p.txt",
                  sep = "\t",
                  row.names = 1,
                  header = T)
view(gt3)
#5p
gt5 <- read.table("~/02_MethodsProject/01_DataFiles/P47_5p/P47_genotype_RHEB_5p.txt",
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

#Add genotyped calls to JLE48 Dataset
#Extract Barcodes
metadata <- onlyJLE48@assays$RNA@counts@Dimnames[[2]]

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
#dfJLE48combined <- left_join(df3p, df5p, by = "barcodes")
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

#Amount of genotyped calls per kit in p47
as.data.frame(table(p47@meta.data$RHEB, p47@meta.data$orig.ident))

####4E####
p47_rheb_combined_tsne <- TSNEPlot(object = p47, split.by = "orig.ident",
                                   label = T, 
                                   cols = c("azure2", "yellow","red", "black"),
                                   pt.size = 2, reduction= "tsene",
                                   group.by = "RHEB", 
                                   label.size = 0.0, 
                                   order = c("alt/alt" ,"alt/ref", "ref/ref","No Call" )) +
  ggtitle("RHEB mutation in Patient 47 GEX") 
p47_rheb_combined_tsne
####3F####
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
