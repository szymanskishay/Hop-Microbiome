#Correlations between CQ values and Diaporthe Abundance
directory <- "path/to/working/directory"
setwd(directory)
ps_hop <- readRDS("hop_pure.rds")
ps_hop %>% tax_glom(taxrank = "Species") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) -> hop_sp_RA     # Sort data frame alphabetically by Genus
ps_hop %>% tax_glom(taxrank = "Species") -> hop_sp

otu_table(hop_sp)

dibundance<-t(otu_table(ps_hop))
dia<-dibundance[,c("OTU_2")]

read.delim("hopmap_a.txt", row.names=1)-> hopmapwithcq

sdh<-hopmapwithcq[row.names(hopmapwithcq) %notin% c("Sample113", "Sample13", "Sample15", "Sample6", "Sample14"),]

#sdh<-sample_data(ps_hop)
synthesis<-merge(sdh, dia, by='row.names')
row.names(synthesis)<-synthesis[,1]
synthesis<-synthesis[,2:length(colnames(synthesis))]
synthesis$shannon<-vegan::diversity(otu_table(ps_hop), MARGIN=2, index="shannon")
synthesis$CQ<-replace(synthesis$CQ, synthesis$CQ == "n/a", 40)
synthesis$CQ<-as.numeric(synthesis$CQ)

synthesis%>%dplyr::filter(Tissue=="Rhizome")->synR
synthesis%>%dplyr::filter(Tissue=="Cone")->synCone

cor.test(x = synCone$OTU_2, y = synCone$CQ) -> glomCone
cor.test(x = synR$OTU_2, y = synR$CQ) -> glomR
cor.test(x = synthesis$OTU_2, y = synthesis$CQ)

summary(lm(CQ~OTU_2*Tissue, synthesis))
