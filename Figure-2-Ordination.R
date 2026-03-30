# Code used to generate ordinations for Figure 2
#108, 109, 118, 119
#### Fungal
directory <- "/path/to/working/directory"
setwd(directory)
ps_hop<-readRDS(file="hop_pure.rds")
sample_data(ps_hop)<-sample_data(read.delim("hopmap.txt", row.names=1))
otu_table(ps_hop) <- prune_taxa(taxa_sums(ps_hop)>10, ps_hop)
ps_hop_placeholder<-ps_hop
ps_hop<-ps_hop_placeholder
otu_table(ps_hop) <-subset(otu_table(ps_hop), select = -c(Sample108, Sample109, Sample118, Sample119)) # problematic

library(metagenomeSeq)
library(ggplot2)
ps_hop.n = phyloseq_to_metagenomeSeq(ps_hop) #n for normalized
p_biom_course<-cumNormStat(ps_hop.n)
biom_quant_course<-cumNorm(ps_hop.n, p=p_biom_course)
normFactors(biom_quant_course)
ps_hop.nf <-MRcounts(biom_quant_course, norm=T)
#create physeq object with normalized otu table
otu_table(ps_hop) <- otu_table(ps_hop.nf, taxa_are_rows = TRUE)
course.ord = ordinate(ps_hop, method ="NMDS", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
course.pca.Tf = plot_ordination(ps_hop, course.ord, color="Tissue") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  #scale_color_manual(values=c("Cone"="lightgreen", "Rhizome"="darkred"))+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1)+
  theme(legend.position = "bottom")+
  labs(tag = "A", title="Tissue")

course.pca.Tf
#####
##### By Week
ps_hop<-readRDS(file="hop_pure.rds")
otu_table(ps_hop) <- prune_taxa(taxa_sums(ps_hop)>10, ps_hop)
ps_cone<-subset_samples(ps_hop, Tissue%in%c("Cone"))
ps_cone.n = phyloseq_to_metagenomeSeq(ps_cone) #n for normalized
p_biom_course<-cumNormStat(ps_cone.n)
biom_quant_course<-cumNorm(ps_cone.n, p=p_biom_course)
normFactors(biom_quant_course)
ps_cone.nf <-MRcounts(biom_quant_course, norm=T)
#create physeq object with normalized otu table
otu_table(ps_cone) <- otu_table(ps_cone.nf, taxa_are_rows = TRUE)
cone.ord = ordinate(ps_cone, method ="NMDS", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
sample_data(ps_cone)$Week<-as.factor(sample_data(ps_cone)$Week)
cone.pca.Wf = plot_ordination(ps_cone, cone.ord, color="Week") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1)+
  theme(legend.position = "bottom")+
  labs(tag = "B", title="Week")

##### Bacterial
directory <- "/path/to/working/directory"
setwd(directory)
ps_hop<-readRDS(file="hop_pure_16.rds")
sample_data(ps_hop)<-sample_data(read.delim("hopmap.txt", row.names=1))
otu_table(ps_hop) <- prune_taxa(taxa_sums(ps_hop)>10, ps_hop)
ps_hop_placeholder<-ps_hop
ps_hop<-ps_hop_placeholder
otu_table(ps_hop) <-subset(otu_table(ps_hop), select = -c(Sample72,)) # problematic
ps_hop <- subset_samples(ps_hop, Week %notin% c("1", "2", "5")) # remove problematic weeks
library(metagenomeSeq)
library(ggplot2)
ps_hop.n = phyloseq_to_metagenomeSeq(ps_hop) #n for normalized
p_biom_course<-cumNormStat(ps_hop.n)
biom_quant_course<-cumNorm(ps_hop.n, p=p_biom_course)
normFactors(biom_quant_course)
ps_hop.nf <-MRcounts(biom_quant_course, norm=T)
#create physeq object with normalized otu table
otu_table(ps_hop) <- otu_table(ps_hop.nf, taxa_are_rows = TRUE)
course.ord = ordinate(ps_hop, method ="NMDS", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
course.pca.Tb = plot_ordination(ps_hop, course.ord, color="Tissue") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  #scale_color_manual(values=c("Cone"="lightgreen", "Rhizome"="darkred"))+
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1)+
  theme(legend.position = "bottom")+
  labs(tag = "C", title="Tissue")

course.pca.T
#####
##### By Week
ps_hop<-readRDS(file="hop_pure_16.rds")
otu_table(ps_hop) <- prune_taxa(taxa_sums(ps_hop)>10, ps_hop)
ps_cone<-subset_samples(ps_hop, Tissue%in%c("Cone"))
ps_cone.n = phyloseq_to_metagenomeSeq(ps_cone) #n for normalized
p_biom_course<-cumNormStat(ps_cone.n)
biom_quant_course<-cumNorm(ps_cone.n, p=p_biom_course)
normFactors(biom_quant_course)
ps_cone.nf <-MRcounts(biom_quant_course, norm=T)
#create physeq object with normalized otu table
otu_table(ps_cone) <- otu_table(ps_cone.nf, taxa_are_rows = TRUE)
cone.ord = ordinate(ps_cone, method ="NMDS", distance="bray")
# this function will visualize the ordination along the two axes that account for the most variation
sample_data(ps_cone)$Week<-as.factor(sample_data(ps_cone)$Week)
cone.pca.Wb = plot_ordination(ps_cone, cone.ord, color="Week") + 
  theme_classic()+
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  theme(plot.title = element_text(size = 18, hjust = 0.5), legend.title=element_text(size=18), legend.text=element_text(size=16)) +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size=16), plot.tag = element_text(size=16, face="bold"))+
  geom_point(size=3, alpha=0.9) + 
  geom_point(size=3, shape=1, color="Black")+
  stat_ellipse(type="norm", linetype=1)+
  theme(legend.position = "bottom")+
  labs(tag = "D", title="Week")
#Figure 2
ggarrange(course.pca.Tf, cone.pca.Wf, course.pca.Tb, cone.pca.Wb)