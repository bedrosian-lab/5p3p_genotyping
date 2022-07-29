#Load combined object
####S3B####
#Get cDNA Length for RTN4 
############### CDNA LENGTH ####
library(tidyverse)
library(GenomicFeatures)
library(scales)

# make 'database' from the gff we used for exome processing
txdb <- makeTxDbFromGFF(".../GRCh38_latest_genomic.gff.gz")

# get transcript lengths for all tx in our txdb
txlens <- transcriptLengths(txdb) %>%
  as.data.frame() %>%
  dplyr::select(tx_name, gene_id, tx_len) %>%
  mutate(tx_gene = paste(tx_name,gene_id, sep = "_")) # make a new column called 'tx_gene' which has unique gene and transcript ID (used for joining later on)

view(txlens)
