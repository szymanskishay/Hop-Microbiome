directory <- "path/to/working/directory"
setwd(directory)
ps_hop<-readRDS(file="hop_pure.rds")
sample_data(ps_hop)<-sample_data(read.delim("hopmap.txt", row.names=1))
otu_table(ps_hop) <- prune_taxa(taxa_sums(ps_hop)>10, ps_hop) # filter low abundance problems


library(metagenomeSeq) # for normalization to help with downstream calculations
ps_hop.n = phyloseq_to_metagenomeSeq(ps_hop) #n for normalized
p_biom_ps_hop<-cumNormStat(ps_hop.n)
biom_quant_ps_hop<-cumNorm(ps_hop.n, p=p_biom_ps_hop)
normFactors(biom_quant_ps_hop)
ps_hop.nf <-MRcounts(biom_quant_ps_hop, norm=T) 
otu_ps_hop <- as.data.frame(otu_table(ps_hop))
taxa_ps_hop <- as.data.frame(as.matrix(tax_table(ps_hop)))
metadata_ps_hop <- as.data.frame(as.matrix(sample_data(ps_hop)))

adonis.whole.T<-adonis2(t(otu_ps_hop) ~ Tissue,  data=metadata_ps_hop, permutations=9999) #Models relationship between Tissue and normalized OTU table
vegan::vegdist(t(otu_ps_hop), method="bray") -> dist_otu_hop #make base info
permdisp_otu_hop.tissue <- betadisper(dist_otu_hop, metadata_ps_hop$Tissue)#Bdisp for tissue
betadisp.spray.tissue<-anova(permdisp_otu_hop.tissue, permutations = 9999)#set to write
######
#reset ps_hop 
ps_hop<-readRDS(file="hop_pure.rds")
sample_data(ps_hop)<-sample_data(read.delim("hopmap.txt", row.names=1))
otu_table(ps_hop) <- prune_taxa(taxa_sums(ps_hop)>10, ps_hop)
##### Get effect by week (restricted to the cones)
ps_cone<-subset_samples(ps_hop, Tissue%in%c("Cone"))
ps_cone.n = phyloseq_to_metagenomeSeq(ps_cone) #n for normalized
p_biom_ps_cone<-cumNormStat(ps_cone.n)
biom_quant_ps_cone<-cumNorm(ps_cone.n, p=p_biom_ps_cone)
normFactors(biom_quant_ps_cone)
ps_cone.nf <-MRcounts(biom_quant_ps_cone, norm=T) 
otu_ps_cone <- as.data.frame(otu_table(ps_cone))
taxa_ps_cone <- as.data.frame(as.matrix(tax_table(ps_cone)))
metadata_ps_cone <- as.data.frame(as.matrix(sample_data(ps_cone)))

adonis.cone.week<-adonis2(t(otu_ps_cone) ~ Week, data=metadata_ps_cone, permutations=9999) 
vegan::vegdist(t(otu_ps_cone), method="bray") -> dist_otu_hop #make base info
permdisp_otu_hop.Week <- betadisper(dist_otu_hop, metadata_ps_cone$Week)#Bdisp for Week
betadisp.spray.Week<-anova(permdisp_otu_hop.Week, permutations = 9999)#set to write

# This gets the Fungal side of Table 1
