# Contact Shay Szymanski for any questions on running this code
# shemans6@msu.edu
library(devtools)
BiocManager::install("DESeq2")
BiocManager::install("genefilter")
install.packages("pheatmap")
BiocManager::install("edgeR")
install.packages("RJSONIO")
install.packages('Tax4Fun_0.3.1.tar.gz', repos=NULL, type='source')
devtools::install_github("joey711/biom")
devtools::install_github("xia-lab/MicrobiomeAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
devtools::install_github("kylebittinger/qiimer")
library(MicrobiomeAnalystR)
# Installing MicrobiomeAnalystR can be difficult, and I cannot assure that this still does it well
directory <- "path/to/directory"
setwd(directory)
ps_hops<-readRDS("hop_pure_16.rds")
'%notin%' <- Negate('%in%')
ps_hops_f <- subset_samples(ps_hops, Week %notin% c("1", "2", "3", "4", "5", "55"))
mbSet<-Init.mbSetObj()
mbSet<-SetModuleType(mbSet, "mdp")
mbSet<-ReadSampleTable(mbSet, "sample_online_cone_only.txt");
mbSet<-Read16STaxaTable(mbSet, "tax_online_cone.csv");
mbSet<-Read16SAbundData(mbSet, "otu_online_cone.csv","text","Others/Not_specific","T");
mbSet<-SanityCheckData(mbSet, "text", disableFilter = FALSE);
#mbSet<-SanityCheckSampleData(mbSet);
#mbSet<-SetMetaAttributes(mbSet, "1")
#mbSet<-PlotLibSizeView(mbSet, "norm_libsizes_0","png");
#mbSet<-SanityCheckData(mbSet, "text","count","true");
#mbSet<-SanityCheckSampleData(mbSet);
#mbSet<-SetMetaAttributes(mbSet, "1")
#mbSet<-PlotLibSizeView(mbSet, "norm_libsizes_1","png");
mbSet<-CreatePhyloseqObj(mbSet, "text","Others/Not_specific","F")
mbSet<-ApplyAbundanceFilter(mbSet, "prevalence", 5, 0.1);
mbSet<-ApplyVarianceFilter(mbSet, "iqr", 0.1);
#mbSet<-GetLibscale(mbSet);
mbSet<-PerformNormalization(mbSet, "none", "none", "none");
mbSet<-PerformLefseAnal(mbSet, 0.1, "fdr", 2.0, "Tissue","F","NA","OTU");
mbSet<-PlotLEfSeSummary(mbSet, 15, "dot", "SupplementalFigure3","png");
#mbSet<-PerformLefseAnal(mbSet, 0.05, "fdr", 2.0, "Week","F","NA","OTU");

# Supplemental Figure 3 was Manually edited output for displaying taxonomic names

#
##### WEEK ####

mbCone<-Init.mbSetObj()
mbCone<-SetModuleType(mbCone, "mdp")
mbCone<-ReadSampleTable(mbCone, "extra_map_cone.txt");
mbCone<-Read16STaxaTable(mbCone, "extra_taxtable_cone.txt");
mbCone<-Read16SAbundData(mbCone, "extra_OTU_table_cone.txt","text","Others/Not_specific","T","false");
mbCone<-SanityCheckData(mbCone, "text","sample","true");
mbCone<-SanityCheckSampleData(mbCone);
mbCone<-SetMetaAttributes(mbCone, "1")
mbCone<-PlotLibSizeView(mbCone, "norm_libsizes_0","png");
mbCone<-SanityCheckData(mbCone, "text","count","true");
mbCone<-SanityCheckSampleData(mbCone);
mbCone<-SetMetaAttributes(mbCone, "1")
mbCone<-PlotLibSizeView(mbCone, "norm_libsizes_1","png");
mbCone<-CreatePhyloseqObj(mbCone, "text","Others/Not_specific","F" , "false")
mbCone<-ApplyAbundanceFilter(mbCone, "prevalence", 5, 0.1);
mbCone<-ApplyVarianceFilter(mbCone, "iqr", 0.1);
mbCone<-PerformNormalization(mbCone, "none", "none", "none", "true",1000)
mbCone<-PerformLefseAnal(mbCone, 0.1, "fdr", 2.0, "Week", "F", "NA", "OTU")

mbSet[["analSet"]][["lefse"]][["resTable"]]->conetable

cone_plot<-conetable[order(conetable$LDAscore),]
cone_plot$Name=factor(row.names(cone_plot), levels=row.names(cone_plot))
cone_plot %>% filter(abs(LDAscore) >=2) %>% filter(FDR < 0.05)-> cone_plot 
cone_plot_tis <- cone_plot
cone_plot_tis %>% mutate(avg=(rowSums(cone_plot_tis[3:4])/2)) %>%
  mutate(ConePor = Cone/rowSums(across(c("Cone", "Rhizome")))) %>%
  mutate(RhizPor = Rhizome/rowSums(across(c("Cone", "Rhizome")))) -> cone_plot_tis2


mean(cone_plot_tis[1,3], cone_plot_tis[1,4])


cone_plot %>% mutate(avg=(rowSums(cone_plot[3:7])/5))%>%
  mutate(stdev = sqrt(((cone_plot[,3]-avg)^2 + (cone_plot[,4]-avg)^2 + (cone_plot[,5]-avg)^2 + (cone_plot[,6]-avg)^2 + (cone_plot[,7]-avg)^2)))%>%
  mutate(Six = (cone_plot[,3]-avg)/stdev) %>%
  mutate(Seven = (cone_plot[,4]-avg)/stdev) %>%
  mutate(Eight = (cone_plot[,5]-avg)/stdev) %>%
  mutate(Nine = (cone_plot[,6]-avg)/stdev) %>%
  mutate(Ten = (cone_plot[,7]-avg)/stdev) -> cone_plot_test

cone_plot

coreplot<-ggplot(cone_plot_test, aes(y=factor(Name), x=LDAscore))+
  theme_classic()+
  geom_point(stat="identity", size=4)+
  labs(y="Species", x="LDAscore")+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major.y = element_line(color = 'grey', linetype=5),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(vjust = 1),
        legend.title=element_blank(),
        legend.margin = margin(0, -10, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
library(tidyverse)
cone_plot_test %>% dplyr::select(Name, Six, Seven, Eight, Nine, Ten) %>%
  pivot_longer(!Name, names_to="Week")-> cone_data_3

cone_plot_tis2 %>% dplyr::select(Name, ConePor, RhizPor) %>%
  pivot_longer(!Name, names_to="Tissue") -> Tissue_data

Tissue_data$Tissue <- gsub("ConePor", "Cone", Tissue_data$Tissue)
Tissue_data$Tissue <- gsub("RhizPor", "Rhizome", Tissue_data$Tissue)
plotside<-ggplot(cone_data_3,aes(x=Week, y=factor(Name), fill=value))+
  theme_classic()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_blank(),  axis.line = element_blank())+
  geom_tile(aes(fill = value, width = 0.9, height=0.9), size=0.5, color="black")+
  scale_fill_gradient2(limits=c(min(cone_data_3$value), max(cone_data_3$value)),
                       breaks=c(min(cone_data_3$value), 0, max(cone_data_3$value)),
                       low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       label = function(x) sprintf("%.2f",x))+
  #scale_x_discrete(position="top", limits=c("Cone", "Rhizome"))+
  scale_x_discrete(position="top", limits=c("Six", "Seven", "Eight", "Nine", "Ten"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=33, hjust=0.2, size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title.align = 0.5,
        plot.background = element_blank(),
        legend.position="right")+
  labs(fill="Deviation from Mean")
library(ggpubr)
tiff(filename="Supplemental-Figure-4.tiff", width = 1000, height = 500, units="px")
ggarrange(coreplot,NULL ,plotside, nrow=1, common.legend = FALSE, legend="right", align='hv', widths = c(1,-0.3,1))
dev.off()

#Supplemental Figure 4 was manually edited to have taxonomic assignments displayed instead of the ZOTU


