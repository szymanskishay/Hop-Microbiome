# Importing and cleaning the data for use in later scripts 
library(phyloseq)
library(Biostrings)
library(data.table)
library(tidyverse)


##### 
options(scipen = 999) # sets the number of decimals displayed on screen
set.seed(0413) # setting the seed for reproducibility so analyses can be reproduced [0413]
directory <- "path/to/working/directory"
setwd(directory)

fungi_otus<- read.delim("otutable_UPARSE_225bp.txt", row.names = 1) #import OTU table from UPARSE
fungi_otus_phy <-otu_table(fungi_otus,
                           taxa_are_rows = TRUE) # OTU table for phyloseq object
fungi_metadata <-read.delim("hopmap.txt",
                            row.names=1) # read metadata

fungi_metadata_phy <-sample_data(fungi_metadata) # metadata for phyloseq object
fungi_seqs<-readDNAStringSet("otus_225bp.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE) # read DNA sequences for phyloseq
fungi_taxonomy<-read.delim("constax_taxonomy.txt",
                           header=TRUE, 
                           row.names=1) # read in taxonomy information to dataframe
fungi_taxonomy_phy <- tax_table(as.matrix(fungi_taxonomy)) # put into phyloseq readable format
hops<-phyloseq(fungi_metadata_phy, fungi_otus_phy, fungi_taxonomy_phy, fungi_seqs) #make into phyloseq object
df.mainseq<-as.data.frame(sample_data(fungi_metadata_phy))# Convert to data frame for exporting

real<-subset_taxa(hops, Kingdom!="Viridiplantae") # keeps all taxa that are not Viridiplantae
real<-subset_taxa(real, Kingdom!="Metazoa") # keeps all taxa that are not Metazoa
hop_pure<-subset_taxa(real, Family!="Malasseziaceae") # keeps all taxa that are not Malasseziaceae (Common skin microflora)

df_hop_pure<-as.data.frame(sample_data(hop_pure)) # converts sample data to a data frame
df_hop_pure$LibrarySize_hop_pure <- sample_sums(hop_pure) #Add a column for Library Size
df_hop_pure <- df_hop_pure[order(df_hop_pure$LibrarySize_hop_pure),] #sorts by library size
df_hop_pure$Index <- seq(nrow(df_hop_pure)) #Adds another layer of order
write.csv(df_hop_pure, file = "rank_sums_hop_pure.csv") # writes a .csv of all samples organized by library size. Used to exclude samples.

df_hops_pure <-as.data.frame(sample_data(hop_pure))

otu_table(hop_pure) <- subset(otu_table(hop_pure), select = -c(Sample113,
                                                               Sample13,
                                                               Sample15,
                                                               Sample6,
                                                               Sample14))
#remove taxa with less than 1000 reads post-plant/metazoa/malasezzia removal

saveRDS(hop_pure, file="hop_pure.rds", compress=FALSE) # creates an R object for easy import to other scripts, avoiding the need to purify each time you start a session.

#Export for each facet of the reduced phyloseq object, if you please. 
write.table(sample_data(hop_pure), "extra_map.txt")
write.table(otu_table(hop_pure), "extra_OTU_table.txt")
write.table(tax_table(hop_pure), "extra_taxtable.txt")
