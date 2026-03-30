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
remotes::install_github("ropensci/taxize")
devtools::install_github("grunwaldlab/metacoder")
devtools::install_github("xia-lab/MicrobiomeAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
devtools::install_github("kylebittinger/qiimer")
library(MicrobiomeAnalystR)
# Installing MicrobiomeAnalystR can be difficult, and I cannot assure that this still does it well
directory <- "path/to/directory"
setwd(directory)
ps_hops<-readRDS("hop_pure.rds")
mbSet<-Init.mbSetObj()
mbSet<-SetModuleType(mbSet, "mdp")
mbSet<-ReadSampleTable(mbSet, "extra_map.txt");
mbSet<-Read16STaxaTable(mbSet, "extra_taxtable.txt");
mbSet<-Read16SAbundData(mbSet, "extra_OTU_table.txt","text","Others/Not_specific","T");
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
mbSet<-PlotLEfSeSummary(mbSet, 15, "dot", "Figure3","png");
#mbSet<-PerformLefseAnal(mbSet, 0.05, "fdr", 2.0, "Week","F","NA","OTU");

## Figure 3 was then manually edited to display the taxonomic assignment of the OTU instead of the OTU itself


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
mbCone<-CreatePhyloseqObj(mbCone, "text","Others/Not_specific","F")
mbCone<-ApplyAbundanceFilter(mbCone, "prevalence", 5, 0.1);
mbCone<-ApplyVarianceFilter(mbCone, "iqr", 0.1);
mbCone<-PerformNormalization(mbCone, "none", "none", "none")
mbCone<-PerformLefseAnal(mbCone, 0.1, "fdr", 2.0, "Week", "F", "NA", "OTU")
mbCone<-PlotLEfSeSummary(mbCone, 15, "dot", "bar_graph_Cone","png");

mbCone[["analSet"]][["lefse"]][["resTable"]]->conetable

cone_plot<-conetable[order(conetable$LDAscore),]
cone_plot$Name=factor(row.names(cone_plot), levels=row.names(cone_plot))
cone_plot %>% filter(LDAscore >=2) %>% relocate('10', .after = '9')-> cone_plot 
cone_plot %>% mutate(avg=(rowSums(cone_plot[3:12])/10))%>%
  mutate(stdev = sqrt(((cone_plot[,3]-avg)^2 + (cone_plot[,4]-avg)^2 + (cone_plot[,5]-avg)^2 + (cone_plot[,6]-avg)^2 + (cone_plot[,7]-avg)^2 + (cone_plot[,8]-avg)^2 + (cone_plot[,9]-avg)^2+ (cone_plot[,10]-avg)^2+ (cone_plot[,11]-avg)^2+ (cone_plot[,12]-avg)^2)/10))%>%
  mutate(One = (cone_plot[,3]-avg)/stdev) %>%
  mutate(Two = (cone_plot[,4]-avg)/stdev) %>%
  mutate(Three = (cone_plot[,5]-avg)/stdev) %>%
  mutate(Four = (cone_plot[,6]-avg)/stdev) %>%
  mutate(Five = (cone_plot[,7]-avg)/stdev) %>%
  mutate(Six = (cone_plot[,8]-avg)/stdev) %>%
  mutate(Seven = (cone_plot[,9]-avg)/stdev) %>%
  mutate(Eight = (cone_plot[,10]-avg)/stdev) %>%
  mutate(Nine = (cone_plot[,11]-avg)/stdev) %>%
  mutate(Ten = (cone_plot[,12]-avg)/stdev) -> cone_plot_test

mbSet[["dataSet"]][["taxa_table"]]->taxtab
taxtab %>% as.data.frame(.) -> taxtab
taxtab$Name<- row.names(taxtab)
CPT<-merge(cone_plot_test, taxtab, by="Name")
CPT %>% filter(LDAscore > 2) %>% filter(FDR < 0.05) -> CPTa

coreplot<-ggplot(CPTa, aes(y=factor(Species), x=LDAscore))+
  theme_classic()+
  geom_point(stat="identity", size=4)+
  labs(y="Species", x="LDAscore")+
  theme(legend.position = "none",
        plot.margin = unit(c(0,0.1,0,0), 'lines'),
        panel.grid.major.y = element_line(color = 'grey', linetype=5),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.title=element_blank(),
        legend.margin = margin(0, -10, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
library(tidyverse)
cone_plot_test %>% dplyr::select(Name, One, Two, Three, Four, Five, Six, Seven, Eight, Nine, Ten) %>%
  pivot_longer(!Name, names_to="Week")-> cone_data_3
plotside<-ggplot(cone_data_3,aes(x=Week, y=factor(Name), fill=value))+
  theme_classic()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_blank(),  axis.line = element_blank())+
  geom_tile(aes(fill = value, width = 0.7, height=0.7), size=0.5, color="black")+
  scale_fill_gradient2(limits=c(min(cone_data_3$value), max(cone_data_3$value)),
                       breaks=c(min(cone_data_3$value), 0, max(cone_data_3$value)),
                       low = "darkblue",
                       mid = "#eeee00",
                       high = "darkred",
                       label = function(x) sprintf("%.2f",x))+
  scale_x_discrete(position="top", limits=c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=33, hjust=0.2, size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title.align = 0.5,
        plot.background = element_blank(),
        legend.position="right")+
  labs(fill=" Deviation \n from Mean")
library(ggpubr)
tiff(filename="Figure_5.tiff", width = 1000, height = 600, units="px")
ggarrange(coreplot, plotside, nrow=1, common.legend = FALSE, legend="right", align='hv', widths = c(2, 2))
dev.off()

# Figure 5 was edited to display the name of the taxonomic assignment instead of OTU
