
# Author: Sophia Pallenberg, Marie-Madlen Pust
# Last updated: 24 February 2022
# Content: Relative and absolute abundances, correlation analyses 

#load libraries
library(readr)
library(readxl)
library(dplyr)
library(vegan)
library(asbio)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(seriation)
library(repr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(metacoder)
library(taxa)
library(ggsignif)
library(magrittr)
library(psych)
library(PerformanceAnalytics)
library(rcompanion)
library(GGally)
library(ggtext)
library(tidyverse)

#set plot options
options(repr.plot.width=12, repr.plot.height = 12)

##Prepare Datafiles and Metadatatable
#load count data and metadata
count_data <- read_excel("count_data_sequins.xlsx", na = c("NA"))
meta_data_04 <- read_excel("meta_data_04.xlsx", na = c("NA"))
#load datafile with rare species <5% of total abundance marked as "others"
count_rareasothers <- read_excel("count_data_sequins_rareasothers.xlsx", sheet = "count_data_sequins_transposed")

# convert to data frame
count_data <- data.frame(count_data)
meta_data_04 <- data.frame(meta_data_04)
count_rareasothers <- as.data.frame(count_rareasothers)

#set metadata columns as numeric:
meta_data_04$time_point <- factor(as.character(meta_data_04$time_point), levels= c("V1", "V2", "V3"))
meta_data_04$bmi <- as.numeric(as.character(meta_data_04$bmi))
meta_data_04$age <- as.numeric(as.character(meta_data_04$age))
meta_data_04$sweat_chloride <- as.numeric(as.character(meta_data_04$sweat_chloride))
meta_data_04$fev1_perc_pred <- as.numeric(as.character(meta_data_04$fev1_perc_pred))
meta_data_04$kaftrio_intake_weeks <- as.numeric(as.character(meta_data_04$kaftrio_intake_week))
meta_data_04$mef25_perc_pred <- as.numeric(as.character(meta_data_04$mef25_perc_pred))
meta_data_04$lci <- as.numeric(as.character(meta_data_04$lci))
meta_data_04$mri_morphology_score <- as.numeric(as.character(meta_data_04$mri_morphology_score))
meta_data_04$mri_global_score <- as.numeric(as.character(meta_data_04$mri_global_score))
meta_data_04$npd_cumdepol <- as.numeric(as.character(meta_data_04$npd_cumdepol))
meta_data_04$npd_mean_sermet <- as.numeric(as.character(meta_data_04$npd_mean_sermet))
meta_data_04$b_adrenergic_sweat_rate <- as.numeric(as.character(meta_data_04$b_adrenergic_sweat_rate))
meta_data_04$delta_bmi_v1 <- as.numeric(as.character(meta_data_04$delta_bmi_v1))
meta_data_04$delta_sweatchloride_v1 <- as.numeric(as.character(meta_data_04$delta_sweatchloride_v1))
meta_data_04$delta_fev1_perc_pred_v1 <- as.numeric(as.character(meta_data_04$delta_fev1_perc_pred_v1))
meta_data_04$delta_mef25_perc_pred_v1 <- as.numeric(as.character(meta_data_04$delta_mef25_perc_pred_v1))
meta_data_04$delta_lci <- as.numeric(as.character(meta_data_04$delta_lci))
meta_data_04$delta_mri_morphology <- as.numeric(as.character(meta_data_04$delta_mri_morphology))
meta_data_04$delta_mri_global <- as.numeric(as.character(meta_data_04$delta_mri_global))
meta_data_04$disease_centile <- as.numeric(as.character(meta_data_04$disease_centile))
meta_data_04$samplecode <- as.numeric(as.character(meta_data_04$samplecode))
meta_data_04$psa_chronic <- as.numeric(as.character(meta_data_04$psa_chronic))

#Set Rownames and Transpose data:
rownames(count_data) <- count_data$Organism
count_data$Organism <- NULL
count_data_t <- data.frame(t(count_data))
rownames(meta_data_04) <- meta_data_04$sample_name

# calculate microbiome parameters and add to meta_data table:
meta_data_04$shannon_div <- vegan::diversity(count_data_t, index = "shannon")
meta_data_04$simpson_div <- vegan::diversity(count_data_t, index = "simpson")
meta_data_04$spec_number <- vegan::specnumber(count_data_t) 
meta_data_04$total_bacterial_load <- rowSums(count_data_t)
meta_data_04$PA_count <- count_data_t$Pseudomonas.aeruginosa
meta_data_04$PA_count_pos <- ifelse(meta_data_04$PA_count > 0, "Pos", "Neg")
meta_data_04$SA_count <- count_data_t$Staphylococcus.aureus
meta_data_04$HiB_count <- count_data_t$Haemophilus.influenzae
meta_data_04$evenness <- meta_data_04$shannon_div / log(meta_data_04$spec_number) 


#Create bargraphs of relative and absolute abundance
rownames(count_rareasothers) <- count_rareasothers$pseudonym
count_rareasothers$pseudonym <- NULL
#add metadata to count data rareasothers table
count_rareasothers$patient_number <- meta_data_04$patient_number
count_rareasothers$sample_type <- meta_data_04$sample_type
count_rareasothers$samplecode <- meta_data_04$samplecode
count_rareasothers$visit <- meta_data_04$time_point
count_rareasothers$genotype <- meta_data_04$mutation_type
count_rareasothers$age <- meta_data_04$age
count_rareasothers$patient_id <- meta_data_04$patient_id
count_rareasothers$PA_count <- count_data_t$Pseudomonas.aeruginosa
count_rareasothers$psa_chronic <- meta_data_04$psa_chronic

#Change character values to numeric codes:
count_rareasothers$visit <- ifelse(meta_data_04$time_point == "V1", 1, ifelse(meta_data_04$time_point == "V2", 2, 3))
count_rareasothers$genotype <- ifelse(meta_data_04$mutation_type == "F/F", 1, ifelse(meta_data_04$mutation_type == "F/MF", 2, NA))
count_rareasothers$sample_type <- ifelse(meta_data_04$sample_type == "throat swab", 1, ifelse(meta_data_04$sample_type == "sputum", 2, NA))

# Remove empty/sterile with no data (RA_Kaf_F17_2005_V2):
count_rareasothers2 <- count_rareasothers[-c(34), ]
rowSums(count_rareasothers2[0:33])

#gather count table to long format
species_long <- gather(count_rareasothers2, key = "species",value = "absolute abundance", 
                       -c(patient_id, patient_number, sample_type, samplecode, patient_number, visit, genotype, age, PA_count, psa_chronic))
species_long$species <- gsub(x = species_long$species, pattern = "\\.", replacement = " ")

#define Visits and Sample Types:
visits <- c("1" = "Induced Sputum V1", "2" = "Throat Swab V1", "3" = "*", "4" = "Throat Swab V2","5" = "Throat Swab V3")
sample <- c("1" = "Throat swab", "2" = "induced sputum")

##Create filled barplot of Species relative Abundance
#define colors for species:
cols <- c("Atopobium parvulum" =	"#67001f",
          "Bifidobacterium longum" = "#9e0142",
          "Capnocytophaga gingivalis"	 ="#d53e4f",
          "Capnocytophaga leadbetteri"  =	"#d53e4f",
          "Capnocytophaga ochracea"	 ="#d53e4f",
          "Capnocytophaga sputigena" =	"#d53e4f",
          "Fusobacterium periodonticum"	 ="#f46d43",
          "Haemophilus influenzae" =	"#fdae61",
          "Haemophilus parainfluenzae" =	"#fdae61",
          "Neisseria subflava" =	"#fee08b",
          "Others" = "#ffffbf",
          "Porphyromonas asaccharolytica"	 ="#e6f598",
          "Porphyromonas gingivalis" =	"#e6f598",
          "Prevotella denticola" =	"#abdda4",
          "Prevotella enoeca" =	"#abdda4",
          "Prevotella fusca" =	"#abdda4",
          "Prevotella jejuni"	 ="#abdda4",
          "Prevotella melaninogenica" =	"#abdda4",
          "Pseudomonas aeruginosa" ="#003c30"	,
          "Rothia aeria" =	"#66c2a5",
          "Rothia dentocariosa" =	"#66c2a5",
          "Rothia mucilaginosa"	 ="#66c2a5",
          "Staphylococcus aureus" = "#762a83",
          "Streptococcus cristatus"	 ="#3288bd",
          "Streptococcus equinus"	 ="#3288bd",
          "Streptococcus mitis" =	"#3288bd",
          "Streptococcus oralis" =	"#3288bd",
          "Streptococcus parasanguinis"	 ="#3288bd",
          "Streptococcus pseudopneumoniae" =	"#3288bd",
          "Streptococcus salivarius" =	"#3288bd",
          "Streptococcus sanguinis" =	"#3288bd",
          "Veillonella atypica" =	"#5e4fa2",
          "Veillonella dispar" =	"#5e4fa2")

#create relative abundance barplot:        
rel_abund <-
  ggplot(species_long, aes(x= reorder (patient_number, age), y=`absolute abundance`, fill=species)) +
  geom_bar(position = "fill", stat = "identity", color = "white", width = 1) +
  facet_grid(~ samplecode, labeller = as_labeller(visits), scales = "free_x", space = "free_x") +
  labs(y="Relative abundance", x="") +
  scale_y_continuous(expand = c(0,0))+
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual("Species", values = cols) +
  theme(legend.position = "right", 
        legend.text =element_text(size=8, face = "italic"), 
        legend.title = element_blank(),
        axis.title = element_text(face = "bold")) 

#create absolute abundance circular plot:      
species_long$absolute_abundance_2 <- log10(species_long$`absolute abundance`+1)
abs_abund <-
  ggplot(species_long, aes(x= reorder (patient_number, age), y=species, colour=species, size=absolute_abundance_2)) +
  geom_point() +
  scale_size(range = c(-1,12)) +
  facet_grid(~ samplecode, labeller = as_labeller(visits), scales = "free_x", space = "free_x") +
  labs(y=" ", x="") +
  theme_light() +
  guides(colour=guide_legend(ncol=1)) +
  scale_colour_manual("Species", values = cols)+
  scale_y_discrete(limits=rev) +
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, size = 8), 
        axis.text.y = element_text(face="italic"),
        plot.title = element_text(size=16, hjust=0.5, face = "bold"), 
        axis.title = element_text(face = "bold")) 

# save absolute abundance plot
pdf("Supplementary_Figure_S01.pdf", width=13, height = 7)
abs_abund
dev.off()

#################################################################################################################
##Create correlation matrices 
# use new meta_data table to format:
meta_data_05 <- meta_data_04

# set characters to numeric values:
meta_data_05$time_point_numb <- ifelse(meta_data_05$time_point == "V1", 1, ifelse(meta_data_05$time_point == "V2", 2, 3))
meta_data_05$sex_numb <- ifelse(meta_data_05$sex == "f", 1, ifelse(meta_data_05$sex == "m", 2, NA))
meta_data_05$mutation_type_numb <- ifelse(meta_data_05$mutation_type == "F/F", 1, ifelse(meta_data_05$mutation_type == "F/MF", 2, NA))
meta_data_05$mutation_code <- ifelse(meta_data_05$mutation_type == "F/F", 1, ifelse(meta_data_05$mutation_type == "F/MF", 2, NA))
meta_data_05$sample_type_numb <- ifelse(meta_data_05$sample_type == "throat swab", 1, ifelse(meta_data_05$sample_type == "sputum", 2, NA))
meta_data_05$pre_kaftrio_modulator_numb <- ifelse(meta_data_05$pre_kaftrio_modulator == "Symkevi", 1, ifelse(meta_data_05$pre_kaftrio_modulator == "Orkambi", 2, 3))

# add species from count data to meta data:
meta_data_05$Rothia.mucilaginosa <- count_data_t$Rothia.mucilaginosa
meta_data_05$Prevotella.melaninogenica <- count_data_t$Prevotella.melaninogenica
meta_data_05$Prevotella.jejuni <- count_data_t$Prevotella.jejuni
meta_data_05$Veillonella.atypica <- count_data_t$Veillonella.atypica
meta_data_05$Prevotella.intermedia <- count_data_t$Prevotella.intermedia
meta_data_05$Veilonelle_parvula <- count_data_t$Veillonella.parvula

# filter for baseline and visit 2:
meta_data_05_filter <- meta_data_05
meta_data_05_filter <- meta_data_05_filter[!meta_data_05_filter$time_point=="V3",]
meta_data_05_filter$mutation_code <- as.numeric(as.character(meta_data_05_filter$mutation_code))
meta_data_05_filter$sex <- NULL
meta_data_05_filter$sex_numb <- NULL

# calculate pairwise correlations:
# select parameters for correlation analysis
meta_data_correlation <- select(meta_data_05_filter, 
                                c("age", "bmi", "disease_centile",
                                  "fev1_perc_pred", "mef25_perc_pred", "lci",
                                  "sweat_chloride","npd_mean_sermet", "b_adrenergic_sweat_rate",
                                  "psa_chronic",  "PA_count", "SA_count", "HiB_count", 
                                  "shannon_div", "simpson_div", "spec_number", "evenness", "total_bacterial_load",   
                                  "Rothia.mucilaginosa", "Prevotella.melaninogenica", 
                                  "Prevotella.jejuni", "Prevotella.intermedia", "Veillonella.atypica", "Veilonelle_parvula",)) 

# create a  Correlation matrix:
meta_data_correlation_matrix <-data.matrix(meta_data_correlation) 
correlation_matrix <- cor(meta_data_correlation_matrix, method="spearman", use="na.or.complete")

# calculate p values
colnames(meta_data_correlation) <- c("Age", "BMI", "Disease centile", "FEV1", "MEF25", 
                                     "LCI2.5", "Sweat chloride concentration", "NPD Sermet score", "b-adrenergic sweat rate", 
                                     "Chronic P.aeruginosa colonization", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Haemophilus influenzae", 
                                     "Shannon diversity index", "Simpson diversity index", "Species number", "Pielou's evenness index", "Total bacterial load", 
                                     "Rothia mucilaginosa", "Prevotella melaninogenica", "Prevotella jejuni", "Prevotella intermedia",
                                     "Veillonella atypica", "Veillonella parvula")
testRes = cor.mtest(meta_data_correlation, method="spearman",conf.level = 0.95)

# set Column and Row Names:
colnames(correlation_matrix) <- c("Age", "BMI", "Disease centile", "FEV1", "MEF25", 
                                  "LCI2.5", "Sweat chloride concentration", "NPD Sermet score", "b-adrenergic sweat rate", 
                                  "Chronic P.aeruginosa colonization", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Haemophilus influenzae", 
                                  "Shannon diversity index", "Simpson diversity index", "Species number", "Pielou's evenness index", "Total bacterial load", 
                                  "Rothia mucilaginosa", "Prevotella melaninogenica", "Prevotella jejuni", "Prevotella intermedia",
                                  "Veillonella atypica", "Veillonella parvula")
rownames(correlation_matrix) <- c("Age", "BMI", "Disease centile", "FEV1", "MEF25", 
                                  "LCI2.5", "Sweat chloride concentration", "NPD Sermet score", "b-adrenergic sweat rate", 
                                  "Chronic P.aeruginosa colonization", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Haemophilus influenzae", 
                                  "Shannon diversity index", "Simpson diversity index", "Species number", "Pielou's evenness index", "Total bacterial load", 
                                  "Rothia mucilaginosa", "Prevotella melaninogenica", "Prevotella jejuni", "Prevotella intermedia",
                                  "Veillonella atypica", "Veillonella parvula")

# create and store correlation plot
pdf(file="Figure_07.pdf")
par(xpd = TRUE) 
corrplot(correlation_matrix, type = 'upper', method = "color", tl.col = 'black', tl.cex = 0.7, tl.srt = 90, #order="hclust", hclust.method = "ward.D2",
         mar = c(2, 0, 0, 0), p.mat = testRes$p, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.7, pch.col = 'black', diag = FALSE,
         win.asp = 1) -> p
dev.off()

# calculate pairwise correlations with ggpairs:
# create data.frame with parameters of interest:
# remove gender from metadata table
meta_data_05$sex <- NULL 
meta_data_05$sex_numb <- NULL
meta_data_cor_2 <- select(meta_data_05, 
                          c("time_point", "age", "bmi", "disease_centile", "mutation_type", "psa_chronic",
                            "fev1_perc_pred", "mef25_perc_pred", "lci",
                            "sweat_chloride","npd_mean_sermet", "b_adrenergic_sweat_rate",
                            "psa_chronic",  "PA_count", "SA_count", "HiB_count", 
                            "shannon_div", "simpson_div", "spec_number", "evenness", "total_bacterial_load",   
                            "Rothia.mucilaginosa", "Prevotella.melaninogenica", "Prevotella.jejuni", 
                            "Prevotella.intermedia", "Veillonella.atypica", "Veilonelle_parvula",))
meta_data_cor_2$time_point <- str_replace(meta_data_cor_2$time_point, "V1", "Baseline")

#Create and store GGpairs plot comparing lung function and microbial parameters:
pdf("Supplementary_Figure_S02.pdf", height=13, width=14)
ggpairs(meta_data_cor_2, column = c("age", "bmi", "fev1_perc_pred", "mef25_perc_pred", "lci", 
                                    "shannon_div", "simpson_div", "spec_number",  "evenness"), 
        aes(color = time_point, alpha = 0.5 ), 
        axisLabels = "show",
        lower = list(continuous = wrap("smooth", alpha =1, se = FALSE, displayGrid = FALSE)),
        upper = list(continuous = wrap("cor", size = 4, display_grid = FALSE, method = "spearman")), 
        title = "Correlation of clinical data, lung function and airway metagenome parameters",
        theme(plot.title = element_text(hjust = 0.5, face= "bold", size=18)),
        columnLabels = c("Age", "BMI", "FEV1", "MEF25", "LCI2.5", "Shannon Diversity", "Simpson Diversity", "Species Richness", "Pielou's evenness"))
dev.off()

#Create and store GGpairs plot comparing CFTR biomarkers and microbial parameters:
pdf("Supplementary_Figure_S03.pdf", height=13, width=13)
ggpairs(meta_data_cor_2, column = c("age",  "sweat_chloride", "npd_mean_sermet", "b_adrenergic_sweat_rate",
                                    "shannon_div", "simpson_div", "spec_number",  "evenness"), 
        aes(color = time_point, alpha = 0.5 ), 
        lower = list(continuous = wrap("smooth", alpha =1, se = FALSE, displayGrid = FALSE)),
        upper = list(continuous = wrap("cor", size = 4, display_grid = FALSE, method = "spearman")), 
        theme(plot.title = element_blank()), 
        columnLabels = c("Age", "Sweat chloride conc.", "NPD Sermet score", "b-adrenergic sweat rate", 
                         "Shannon Diversity", "Simpson Diversity", "Species Richness", "Pielou's evenness"))
dev.off()


