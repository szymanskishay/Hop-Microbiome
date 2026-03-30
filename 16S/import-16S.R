# Importing and cleaning the data for use in later scripts 
library(phyloseq)
library(Biostrings)
library(data.table)
library(tidyverse)

##### 
options(scipen = 999) # sets the number of decimals displayed on screen 0
set.seed(0413) # setting the seed for reproducibility so analyses can be reproduced [0413]
directory <- "path/to/working/directory"
setwd(directory)

bac_otus<- read.delim("otutable_UNOISE_240bp.txt", row.names = 1) # import zOTU table from UNOISE
bac_otus_phy <-otu_table(bac_otus,
                           taxa_are_rows = TRUE) # read zOTU table for phyloseq
bac_metadata <-read.delim("hopmap_16.txt",
                            row.names=1) # read metadata

bac_metadata_phy <-sample_data(bac_metadata) # read for Phyloseq
bac_seqs<-readDNAStringSet("asv_240bp.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE) #asv sequences reading
bac_taxonomy<-read.delim("constax_taxonomy.txt",
                           header=TRUE, 
                           row.names=1) # read in taxonomy information to dataframe
bac_taxonomy_phy <- tax_table(as.matrix(bac_taxonomy)) # put into phyloseq readable format
hops<-phyloseq(bac_metadata_phy, bac_otus_phy, bac_taxonomy_phy, bac_seqs)
df.mainseq<-as.data.frame(sample_data(bac_metadata_phy))# Convert to data frame for exporting

real<-subset_taxa(hops, Rank_4!="Chloroplast_1") # remove any sequences assigned as chloroplasts
real<-subset_taxa(real, Rank_5!="Mitochondria_1") # remove any sequences assigned as mitochondria
hop_pure<-real 
df_hop_pure<-as.data.frame(sample_data(hop_pure)) # data frame from new set
df_hop_pure$LibrarySize_hop_pure <- sample_sums(hop_pure) #Add a column for Library Size
df_hop_pure <- df_hop_pure[order(df_hop_pure$LibrarySize_hop_pure),] #sorts
df_hop_pure$Index <- seq(nrow(df_hop_pure)) #Adds another layer of order
write.csv(df_hop_pure, file = "rank_sums_hop_pure.csv") #write csv to get samples ordered by library size after filtering chloroplasts/mitochondria


df_hops_pure <-as.data.frame(sample_data(hop_pure))

otu_table(hop_pure) <- subset(otu_table(hop_pure), select = -c(Sample16,
                                                               Sample9,
                                                               Sample11,
                                                               Sample36,
                                                               Sample2,
                                                               Sample39,
                                                               Sample4,
                                                               Sample3,
                                                               Sample26,
                                                               Sample12,
                                                               Sample27,
                                                               Sample15,
                                                               Sample10,
                                                               Sample32,
                                                               Sample30,
                                                               Sample6,
                                                               Sample8,
                                                               Sample22,
                                                               Sample13,
                                                               Sample23,
                                                               Sample41,
                                                               Sample45,
                                                               Sample47,
                                                               Sample42,
                                                               Sample50,
                                                               Sample44,
                                                               Sample33,
                                                               Sample35,
                                                               Sample37,
                                                               Sample20,
                                                               Sample14,
                                                               Sample29,
                                                               Sample1,
                                                               Sample7,
                                                               Sample28,
                                                               Sample18,
                                                               Sample21,
                                                               Sample40,
                                                               Sample24,
                                                               Sample43,
                                                               Sample31,
                                                               Sample48,
                                                               Sample34,
                                                               Sample38,
                                                               Sample19,
                                                               Sample25,
                                                               Sample106,
                                                               Sample60,
                                                               Sample61,
                                                               Sample67,
                                                               Sample52,
                                                               Sample69,
                                                               Sample51))
#remove taxa with less than 1000 reads post-michondrial and chloroplastic removal
  # Note, this is a lot of samples.
saveRDS(hop_pure, file="hop_pure_16.rds", compress=FALSE) #saves phyloseq object for easier starting point on other scripts

#Export for each facet of the reduced phyloseq object, if you please. 
write.table(sample_data(hop_pure), "extra_map.txt")
write.table(otu_table(hop_pure), "extra_OTU_table.txt")
write.table(tax_table(hop_pure), "extra_taxtable.txt")

