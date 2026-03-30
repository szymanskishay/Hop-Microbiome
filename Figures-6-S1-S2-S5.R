library(phyloseq)
library(dplyr)
library(data.table)
library(ggplot2)

ps_hop<-readRDS(file="hop_pure.rds")

# -------- Filter low counts --------
otu_table(ps_hop) <- otu_table(ps_hop)[which(rowSums(otu_table(ps_hop)) >= 10),]

# -------- Process data --------
ps_hop.bar<- ps_hop %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  arrange(Genus)

dat_ps_hop.bar <- data.table(ps_hop.bar)


# --------Collapse low abundance taxa--------

dat_ps_hop.bar[(Abundance <= 0.05), Genus:= "Other"] 

write.csv(dat_ps_hop.bar, "allgenus.csv")

# -------- Split into the two tissue types --------
tissues <- unique(dat_ps_hop.bar$Tissue)

Cone_abun <- dat_ps_hop.bar[Tissue == tissues[1]]
rhizome_abun <- dat_ps_hop.bar[Tissue == tissues[2]]

# -------- Ordercone order by Week then Rep--------
Cone_abun$Week <- factor(Cone_abun$Week, levels = 1:10)
Cone_abun$Rep  <- factor(Cone_abun$Rep,  levels = 1:11)


# --------Create correct Sample order based on Week + Rep--------
sample_order_dat1 <- Cone_abun %>%
  distinct(Sample, Week, Rep) %>%
  arrange(Week, Rep) %>%
  pull(Sample)

Cone_abun$Sample <- factor(Cone_abun$Sample, levels = sample_order_dat1)

# --------order rhizome by Rep only--------

rhizome_abun$Rep <- factor(rhizome_abun$Rep, levels = 1:15)

sample_order_dat2 <- rhizome_abun %>%
  distinct(Sample, Rep) %>%
  arrange(Rep) %>%
  pull(Sample)

rhizome_abun$Sample <- factor(rhizome_abun$Sample, levels = sample_order_dat2)

# --------Colors for fungal abudance--------

my_colors <- c(
  "Aureobasidium"="#ffff00",
  "Alternaria"="#008282",
  "Codinaea"="#18e700",
  "Fusarium"="#ffd1d1",
  "Cladosporium"="#f0d1ff",
  "Cystofilobasidium"="#77003c",
  "Diaporthe"="#004182",
  "Epicoccum"="#f2a400",	
  "Filobasidium"="#ff4400",	
  "Entoloma"="#a1a100",	
  "Bullera"="#77c350",	
  "Papiliotrema"="#0086EB",
  "Curvularia"="#00d5f2",
  "Dactylonectria"="#5b4f3d",
  "Sporobolomyces"="#a954ff",	
  "Flagelloscypha"="#6a006a",
  "Vishniacozyma"="#416600",	
  "Chaetasbolisia"="#ff0000",
  "Diozegia"="#d29E92",
  "Cercospora"="#6aa6ef",
  "Ganoderma"="#ff00f9",
  "Leptospora"="#fff6d1",
  "Microdochium"="#078446",
  "Mycena"="#CCC3B4",
  "Neoascochyta"="#bb0000",
  "Penicillium"="#9c4dad",
  "Periconia"="#d1ffe3",
  "Plectosphaerella"="#99004d",
  "Pseudopithomyces"="#e800a4",
  "Rhinosporidium"="#d1f9ff",
  "Saccharomyces"="#929292",
  "Sampaiozyma"="#4ac925",
  "Serendipita"="#ffe6d1",
  "Sphaerobolus"="#391e71",
  "Sporidesmiella"="#ff6ff2",
  "Thyridium"="#a25203",
  "Tilletiopsis"="#c39797",
  "Trechispora"="#b7fffd",
  "Trichopeziza"="#df9493",
  "Verticillium"="#f86659",
  "Other"="#d1dcff"
)

# -------- variable to lines to weekly plot --------
group_sizes <- c(9, 7, 10, 10, 10, 10, 10, 11, 10)


# -------- Weekly cone samples (Figure 6) --------
ps_hop.barplot1 = ggplot(Cone_abun, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", color = NA, width = 1) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle(paste("PLP -", tissues[1]))+
  theme_classic()+
  theme(axis.text.x= element_blank())+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 5%) \n") +
  xlab("")+
  geom_vline(xintercept = line_positions + 0.5, color = "black", size = 0.5)

ps_hop.barplot1




# -------- Rhizome abundance (Figure S2) --------
ps_hop.barplot2 = ggplot(rhizome_abun, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", color = NA, width = 1) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle(paste("PLP -", tissues[2]))+
  theme_classic()+
  theme(axis.text.x= element_blank())+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 5%) \n") +
  xlab("")

ps_hop.barplot2





# -------- Bacteria abundance --------

# -------- read in data --------
ps_hop<-readRDS(file="hop_pure_16.rds")


taxBeta<-read.delim("constax_taxonomy_BETA.txt", header=TRUE, row.names=1)
tax_table(ps_hop)<-tax_table(as.matrix(taxBeta))
otu_table(ps_hop) <- otu_table(ps_hop)[which(rowSums(otu_table(ps_hop)) >= 10),]### PCR Errors



# -------- Process data --------
ps_hop.bar_16s <- ps_hop %>%
  tax_glom(taxrank = "Rank_3") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Rank_3)           # Sort data frame alphabetically by Genus

dat_ps_hop.bar_16s <- data.table(ps_hop.bar_16s)


#  -------- Collapse low abundance taxa  --------
dat_ps_hop.bar_16s[(Abundance <= 0.01), Rank_3:= "Other"] 



# -------- SPLIT INTO TWO TISSUES --------
tissues_16s <- unique(dat_ps_hop.bar_16s$Tissue)

Cone_abun_16s <- dat_ps_hop.bar_16s[Tissue == tissues_16s[2]]
rhizome_abun_16s <- dat_ps_hop.bar_16s[Tissue == tissues_16s[1]]
Cone_abun_16s <- Cone_abun_16s[!Cone_abun_16s$Week %in% c(1, 2), ]
# -------- Cone: order by Week then Rep --------

Cone_abun_16s$Week <- factor(Cone_abun_16s$Week, levels = c(5, 6, 7, 8, 9, 10))
Cone_abun_16s$Rep  <- factor(Cone_abun_16s$Rep,  levels = 1:11)

# -------- Create correct Sample order based on Week + Rep--------
sample_order_dat1_16s <- Cone_abun_16s %>%
  distinct(Sample, Week, Rep) %>%
  arrange(Week, Rep) %>%
  pull(Sample)

Cone_abun_16s$Sample <- factor(Cone_abun_16s$Sample, levels = sample_order_dat1_16s)
# -------- Rhizome: order by Rep only--------

rhizome_abun_16s$Rep <- factor(rhizome_abun_16s$Rep, levels = 1:15)

sample_order_dat2_16s <- rhizome_abun %>%
  distinct(Sample, Rep) %>%
  arrange(Rep) %>%
  pull(Sample)

rhizome_abun_16s$Sample <- factor(rhizome_abun_16s$Sample, levels = sample_order_dat2_16s)






#  -------- Rank 3 color  --------
Rank_3_colors <-c("Gammaproteobacteria_1"="#ff4200",
  "Betaproteobacteria_1"="#000056",
  "Saccharimonadia_1" = "#10E0FF",
  "Acidimicrobiia_1"="#4ac925",
  "Acidobacteriae_1"="#00d5f2",
  "Actinobacteria_1"="#f191ff",
  "Alphaproteobacteria_1"="#77003c",
  "Bacilli_1"="#008282",
  "Bacteroidia_1"="#f2a400",
  "Blastocatellia_1"="#ffff00",
  "Kd4-96_1"="#2d5817",
  "Nitrososphaeria_1"="#8a5520",
  "Polyangia_1"="#ffd1d1",
  "Thermoleophilia_1"="#004182",
  "Verrucomicrobiae_1"="#900fff",
  "Viciniamibacteria_1"="#d1dcff",
  "Other"="black",
  "Vicinamibacteria_1" = "#ff7b03")






#  -------- Rank 3 plot (Figure S5) --------
ps_hop.barplot= ggplot(dat_ps_hop.bar_16s, aes(x = Sample, y = Abundance, fill = Rank_3)) + 
  facet_wrap(~Tissue, strip.position = "bottom", scales="free_x") +
  geom_bar(stat = "identity", color = NA, width = 1) +
  scale_fill_manual(values = Rank_3_colors)+
  scale_y_continuous(expand = c(0,0)) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("Hop")+
  theme_classic()+
  theme(axis.text.x= element_blank())+
  #theme(axis.text.x = element_text(size = 0, angle = 0, vjust = 0, hjust = 0)) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 14, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Rank_3 > 1%) \n") +
  xlab("")

plot(ps_hop.barplot) 


Rank_2_color <- c("Abditibacteriota_1" = "#9c4dad",
                  "Acidobacteriota_1"= "#a15000",
                  "Actinobacteriota_1"="#a1a100",
                  "Armatimonadota_1" ="#ff0000",
                  "Bacteroidota_1"="#00d5f2",
                  "Bdellovibrionota_1" = "#38f43d",
                  "Chloroflexi_1"="#f0d1ff",
                  "Crenarcheota_1"="#6a006a",
                  "Cyanobacteria_1"= "#f2a400",
                  "Fcpu426_1" = "#ffd1d1",
                  "Firmicutes_1"="#078446",
                  "Gemmatimonadota_1" = "#ffff01",
                  "Latescibacterota_1" = "#a954ff",
                  "Myxococcota_1"="black",
                  "Nitrospirota_1" = "#d1f9ff",
                  "Patescibacteria_1"="#008282",
                  "Planctomycetota_1"="#77003c",
                  "Proteobacteria_1"="#005682",
                  "Pseudomonadaceae_1" = "#4ce24e",
                  "SAR324_clade(Marine_group_B)_1" = "#bd1864",
                  "Spirochaetota_1" = "#ff2106",
                  "Verrucomicrobiota_1"="#ff6ff2",
                  "Other"="white",
                  "WPS-2_1" = "#000056"
)


ps_hop.barplot.rank2= ggplot(dat_ps_hop.bar_16s, aes(x = Sample, y = Abundance, fill = Rank_2)) + 
  facet_wrap(~Tissue, strip.position = "bottom", scales="free_x") +
  #theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity", color = NA, width = 1) +
  scale_fill_manual(values = Rank_2_color)+
  scale_y_continuous(expand = c(0,0)) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("Hop")+
  theme_classic()+
  theme(axis.text.x= element_blank())+
  #theme(axis.text.x = element_text(size = 0, angle = 0, vjust = 0, hjust = 0)) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 14, face = "bold")) +
  theme(legend.position="right")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Rank_2 > 1%) \n") +
  xlab("")

ps_hop.barplot.rank2


#  -------- funguild analysis --------

#in python
 
#python Guilds_v1.0.py -otu /Users/username/Documents/project/otu_table.txt -db fungi


#  -------- funguild plot (Figure S1) --------

fungal_otu_no_hop_withtax_guilds_matched <- read_excel("fungal_otu_no_hop_withtax.guilds_matched.xlsx", 
                                                       sheet = "data._for_r")
View(fungal_otu_no_hop_withtax_guilds_matched)

guild_color <- c('Pathotroph' =  "#e00707",
                  "Pathotroph-Saprotroph" = "#626262",
                  "Symbiotroph" = "#008282",
                  "Pathotroph-Saprotroph-Symbiotroph" = "#b536da",
                  "Pathotroph-Symbiotroph" = "#008141",
                  "Saprotroph" = "#005682",
                  "Saprotroph-Symbiotroph" = "#a15000")
                  
                  

ggplot(fungal_otu_no_hop_withtax_guilds_matched, 
       aes(fill = `ecological guild`, y = per, x = num)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = guild_color) + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Percent of OTUs (1%) \n") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




