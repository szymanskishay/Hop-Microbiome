# Script for Alpha Diversity calculations and plotting for Figure 1
require(phyloseq)
require(vegan)
require(agricolae)
require(data.table)
require(ggplot2)
require(ggrepel)
require(ggpubr)
directory <- "path/to/working/directory"
setwd(directory)
ps_hop<-readRDS(file="hop_pure.rds")
sample_data(ps_hop)<-sample_data(read.delim("hopmap.txt",row.names=1))

sd<-data.frame(sample_data(ps_hop)) #Retrieving metadata from phyloseq
sd$shannon<-vegan::diversity(otu_table(ps_hop), MARGIN=2, index="shannon") #shannon index for whole data set
sd$spec<-specnumber(otu_table(ps_hop), MARGIN=2) 
sd$Pielou<-(sd$shannon/log(sd$spec)) # Pielou's evenness 


shannon_anova_Tissue<-aov(shannon~Tissue, sd) #ANOVA for relationship between Tissue Type and shannon's index
pielou_anova_tissue<-aov(Pielou~Tissue, sd) #ANOVA for relationship between Tissue Type and Pielou's Evenness

tukey_Shannon_tissue<-HSD.test(shannon_anova_Tissue,  "Tissue", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on tissue type
tukey_pielou_tissue<-HSD.test(pielou_anova_tissue, "Tissue", group=TRUE, unbalanced = TRUE) #group means by grop, unbalanced is to deal with unequal sample sizes in comparisons.

#make into new objects for easier plotting
tk.Shannon.tissue<-cbind(tukey_Shannon_tissue$means[,2:9], tukey_Shannon_tissue$groups[order(row.names((tukey_Shannon_tissue$groups))),])

tk.pielou.tissue<-cbind(tukey_pielou_tissue$means[,2:9], tukey_pielou_tissue$groups[order(row.names((tukey_pielou_tissue$groups))),])
##### Plotting for Alpha-Diversity measurements. 
Shannon.tissue.plot.f<-ggplot(tk.Shannon.tissue, aes(x=row.names(tk.Shannon.tissue),
                                                   ymin=Min,
                                                   lower=Q25,
                                                   middle=Q50,
                                                   upper=Q75,
                                                   ymax=Max,
                                                   fill=row.names(tk.Shannon.tissue)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Shannon's Diversity")+
  labs(tag="A")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Cone"="forestgreen", "Rhizome"="brown"))+
  guides(fill=guide_legend(title="Tissue"))

pielou.tissue.plot.f<-ggplot(tk.pielou.tissue, aes(x=row.names(tk.pielou.tissue),
                                                 ymin=Min,
                                                 lower=Q25,
                                                 middle=Q50,
                                                 upper=Q75,
                                                 ymax=Max,
                                                 fill=row.names(tk.pielou.tissue)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Pielou's Evenness")+
  labs(tag="B")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Cone"="forestgreen", "Rhizome"="brown"))+
  guides(fill=guide_legend(title="Tissue"))



#################
### Plots isolated by tissue type
ps_cone<-subset_samples(ps_hop, Tissue%in%c("Cone"))

sd_C<-data.frame(sample_data(ps_cone)) #Retrieving metadata from phyloseq
sd_C$shannon<-vegan::diversity(otu_table(ps_cone), MARGIN=2, index="shannon") #shannon index for whole data set
sd_C$spec<-specnumber(otu_table(ps_cone), MARGIN=2)
sd_C$Pielou<-(sd_C$shannon/log(sd_C$spec))



shannon_anova_Week<-aov(shannon~Week, sd_C) #asessing variance and the like
pielou_anova_Week<-aov(Pielou~Week, sd_C)

tukey_Shannon_Week<-HSD.test(shannon_anova_Week,  "Week", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on Week type
tukey_pielou_Week<-HSD.test(pielou_anova_Week, "Week", group=TRUE, unbalanced = TRUE)
library(stringr) #need this for string sorting to make the weeks actually be in proper order

#make into new objects for easier plotting
tk.Shannon.Week<-cbind(tukey_Shannon_Week$means[,2:9], tukey_Shannon_Week$groups[str_sort(row.names((tukey_Shannon_Week$groups)),numeric=TRUE),])
tk.pielou.Week<-cbind(tukey_pielou_Week$means[,2:9], tukey_pielou_Week$groups[str_sort(row.names((tukey_pielou_Week$groups)), numeric=TRUE),])

real_weeks<-str_sort(row.names(tk.Shannon.Week), numeric=TRUE)
#### Plots by Week (FIG 4)

Shannon.Week.plot.f<-ggplot(tk.Shannon.Week, aes(x=row.names(tk.Shannon.Week),
                                                   ymin=Min,
                                                   lower=Q25,
                                                   middle=Q50,
                                                   upper=Q75,
                                                   ymax=Max,
                                                   fill=row.names(tk.Shannon.Week)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Week")+
  ylab("Shannon's Diversity")+
  labs(tag="A")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_weeks)+
  scale_x_discrete(limits=real_weeks)+
  guides(fill=guide_legend(title="Week"))
pielou.Week.plot.f<-ggplot(tk.pielou.Week, aes(x=row.names(tk.pielou.Week),
                                                 ymin=Min,
                                                 lower=Q25,
                                                 middle=Q50,
                                                 upper=Q75,
                                                 ymax=Max,
                                                 fill=row.names(tk.pielou.Week)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Week")+
  ylab("Pielou's Evenness")+
  labs(tag="B")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_weeks)+
  scale_x_discrete(limits=real_weeks)+
  guides(fill=guide_legend(title="Week"))
#################

# Repeat of above, was originally organized as two separate scripts, hence redundancy.
# Script for Alpha Diversity calculations and plotting for the 16S zOTUs
directory <- "path/to/working/directory"
setwd(directory)
ps_hop<-readRDS(file="hop_pure_16.rds")
sample_data(ps_hop)<-sample_data(read.delim("hopmap_16.txt",row.names=1))
'%notin%'<- Negate('%in%')#useful way to not have to denote as many entries

ps_hop<-subset_samples(ps_hop,Week%notin%c("1", "2", "3", "4", "5"))

sd<-data.frame(sample_data(ps_hop)) #Retrieving metadata from phyloseq
sd$shannon<-vegan::diversity(otu_table(ps_hop), MARGIN=2, index="shannon") #shannon index for whole data set
sd$spec<-specnumber(otu_table(ps_hop), MARGIN=2) 
sd$Pielou<-(sd$shannon/log(sd$spec)) # Pielou's evenness 


shannon_anova_Tissue<-aov(shannon~Tissue, sd) #ANOVA for relationship between Tissue Type and shannon's index
pielou_anova_tissue<-aov(Pielou~Tissue, sd) #ANOVA for relationship between Tissue Type and Pielou's Evenness

tukey_Shannon_tissue<-HSD.test(shannon_anova_Tissue,  "Tissue", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on tissue type
tukey_pielou_tissue<-HSD.test(pielou_anova_tissue, "Tissue", group=TRUE, unbalanced = TRUE) #group means by grop, unbalanced is to deal with unequal sample sizes in comparisons.

#make into new objects for easier plotting
tk.Shannon.tissue<-cbind(tukey_Shannon_tissue$means[,2:9], tukey_Shannon_tissue$groups[order(row.names((tukey_Shannon_tissue$groups))),])

tk.pielou.tissue<-cbind(tukey_pielou_tissue$means[,2:9], tukey_pielou_tissue$groups[order(row.names((tukey_pielou_tissue$groups))),])
##### Plotting for Alpha-Diversity measurements. 
Shannon.tissue.plot.b<-ggplot(tk.Shannon.tissue, aes(x=row.names(tk.Shannon.tissue),
                                                   ymin=Min,
                                                   lower=Q25,
                                                   middle=Q50,
                                                   upper=Q75,
                                                   ymax=Max,
                                                   fill=row.names(tk.Shannon.tissue)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Shannon's Diversity")+
  labs(tag="C")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Cone"="forestgreen", "Rhizome"="brown"))+
  guides(fill=guide_legend(title="Tissue"))

pielou.tissue.plot.b<-ggplot(tk.pielou.tissue, aes(x=row.names(tk.pielou.tissue),
                                                 ymin=Min,
                                                 lower=Q25,
                                                 middle=Q50,
                                                 upper=Q75,
                                                 ymax=Max,
                                                 fill=row.names(tk.pielou.tissue)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Tissue")+
  ylab("Pielou's Evenness")+
  labs(tag="D")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_manual(values=c("Cone"="forestgreen", "Rhizome"="brown"))+
  guides(fill=guide_legend(title="Tissue"))



#################
### Plots isolated by tissue type
'%notin%'<- Negate('%in%') #useful way to not have to denote as many entries

ps_cone<-subset_samples(ps_hop, Tissue%in%c("Cone") & Week%notin%c("1", "2", "3", "4", "5"))

sd_C<-data.frame(sample_data(ps_cone)) #Retrieving metadata from phyloseq
sd_C$shannon<-vegan::diversity(otu_table(ps_cone), MARGIN=2, index="shannon") #shannon index for whole data set
sd_C$spec<-specnumber(otu_table(ps_cone), MARGIN=2)
sd_C$Pielou<-(sd_C$shannon/log(sd_C$spec))



shannon_anova_Week<-aov(shannon~Week, sd_C) #asessing variance and the like
pielou_anova_Week<-aov(Pielou~Week, sd_C)

tukey_Shannon_Week<-HSD.test(shannon_anova_Week,  "Week", group=TRUE, unbalanced=TRUE) #separation of alpha diversity based on Week type
tukey_pielou_Week<-HSD.test(pielou_anova_Week, "Week", group=TRUE, unbalanced = TRUE)
library(stringr) #need this for string sorting to make the weeks actually be in proper order

#make into new objects for easier plotting
tk.Shannon.Week<-cbind(tukey_Shannon_Week$means[,2:9], tukey_Shannon_Week$groups[str_sort(row.names((tukey_Shannon_Week$groups)),numeric=TRUE),])
tk.pielou.Week<-cbind(tukey_pielou_Week$means[,2:9], tukey_pielou_Week$groups[str_sort(row.names((tukey_pielou_Week$groups)), numeric=TRUE),])

real_weeks<-str_sort(row.names(tk.Shannon.Week), numeric=TRUE)
#### Plots by Week (FIG 4)

Shannon.Week.plot.b<-ggplot(tk.Shannon.Week, aes(x=row.names(tk.Shannon.Week),
                                               ymin=Min,
                                               lower=Q25,
                                               middle=Q50,
                                               upper=Q75,
                                               ymax=Max,
                                               fill=row.names(tk.Shannon.Week)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Week")+
  ylab("Shannon's Diversity")+
  labs(tag="C")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_weeks)+
  scale_x_discrete(limits=real_weeks)+
  guides(fill=guide_legend(title="Week"))
pielou.Week.plot.b<-ggplot(tk.pielou.Week, aes(x=row.names(tk.pielou.Week),
                                             ymin=Min,
                                             lower=Q25,
                                             middle=Q50,
                                             upper=Q75,
                                             ymax=Max,
                                             fill=row.names(tk.pielou.Week)))+
  theme_classic()+
  geom_errorbar(aes(ymin=Min, ymax=Max, y=Q50, width=0.5))+
  geom_boxplot(stat="identity")+
  geom_text(aes(y=Max, label=groups), vjust=-0.25)+
  xlab("Week")+
  ylab("Pielou's Evenness")+
  labs(tag="D")+
  theme(plot.tag=element_text(size=18))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), axis.title.x=element_text(size=14, face="bold"),axis.title.y=element_text(size=14, face="bold"),legend.title=element_text(size=14, face="bold"),  legend.text=element_text(size=12))+
  scale_fill_discrete(limits=real_weeks)+
  scale_x_discrete(limits=real_weeks)+
  guides(fill=guide_legend(title="Week"))

#Figure 1
ggarrange(Shannon.tissue.plot.f, pielou.tissue.plot.f, Shannon.tissue.plot.b, pielou.tissue.plot.b, nrow=2, ncol=2, common.legend = TRUE, legend = "bottom")
#Figure 4
ggarrange(Shannon.week.plot.f, pielou.week.plot.f, Shannon.week.plot.b, pielou.week.plot.b, nrow=2, ncol=2, common.legend = TRUE, legend = "bottom")
