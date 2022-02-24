# Author: Marie-Madlen Pust
# Last updated: 24 February 2022
# Content: Sequins normalisation, random forest, network construction, network analysis

# empty global environment
rm(list = ls())

# list of required packages
req_packages <- c('readr', 'string', 'purrr', 'vegan','ggrepel', 'viridis', 'tidyr', 'Hmisc', 'dplyr', 
                  'igraph', 'matrixStats', 'ggpubr', 'ggthemes', 'hrbrthemes', 'plyr', 'readxl', 
                  'randomForest', 'Boruta', 'stringr')

# function for installing/importing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

# define global variables
sig_level = 0.01 # significant level for network correlation analysis
weight_val_pos = 0.2 
weight_val_neg = -0.2
random_seeds <- sample(1:200, 100, replace = FALSE) # repeat random forest with 100 different seeds set

# load or install R packages
ipak(req_packages)

####################################################################################
# import RPMM normalised data tables
count_data_sequins <- read_delim("count_data_SP_sequins_2022_02_08.csv", delim = ";", escape_double = TRUE, na = "NA", trim_ws = TRUE)
count_data_sequins <- data.frame(count_data_sequins)
count_data_sequins_blank <- count_data_sequins[,grepl("Blank", colnames(count_data_sequins))]
count_data_sequins <- count_data_sequins[,!grepl("Blank", colnames(count_data_sequins))]
count_data_sequins <- count_data_sequins[,!grepl("_T", colnames(count_data_sequins))]

# clean species names
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "chromosome", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "genome", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "complete", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "sequence", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "BAC", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "strain", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "subsp", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "ATCC", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "DSM", "")
count_data_sequins$organism <- str_replace_all(count_data_sequins$organism, "____", "")

# set species names as row names
rownames(count_data_sequins) <- count_data_sequins$organism
count_data_sequins$organism <- NULL

# normalise RPMM values based on artificial spike-ins
count_data_sequins_t <- data.frame(t(count_data_sequins))
count_data_sequins_t$total_count <- rowSums(count_data_sequins_t)
count_data_sequins_t_rel <- count_data_sequins_t[1:2781] / count_data_sequins_t$total_count
count_data_sequins_t_rel$total_count <- NULL
count_data_sequins_norm <- count_data_sequins_t_rel[1:2780] / count_data_sequins_t_rel$Sequins
count_data_sequins_norm$Sequins <- NULL
count_data_sequins_norm_00 <- as.data.frame(sapply( count_data_sequins_norm, as.numeric ))
count_data_sequins_norm_00[count_data_sequins_norm_00 == 0] <- NA
rownames(count_data_sequins_norm_00) <- rownames(count_data_sequins_norm)
count_data_sequins_norm_01 <- data.frame(t(count_data_sequins_norm_00))

# filter species 
count_data_sequins_norm_02 <- count_data_sequins_norm_01[which(rowMeans(!is.na(count_data_sequins_norm_01)) > 0.1), ]
count_data_sequins_norm_02[is.na(count_data_sequins_norm_02)] <- 0
count_data_sequins_norm_03 <- count_data_sequins_norm_02 * 1000 
list_rownames <- rownames(count_data_sequins_norm_03)
count_data_sequins_norm_03$names1 <- unlist(map(str_split(list_rownames, "_"),4))
count_data_sequins_norm_03$names2 <- unlist(map(str_split(list_rownames, "_"),3))
count_data_sequins_norm_03$Genus <- ifelse(count_data_sequins_norm_03$names2 == "1", count_data_sequins_norm_03$names1,
                                           ifelse(count_data_sequins_norm_03$names2 == "2", count_data_sequins_norm_03$names1,
                                                  ifelse(count_data_sequins_norm_03$names2 == "3", count_data_sequins_norm_03$names1, 
                                                         count_data_sequins_norm_03$names2)))
count_data_sequins_norm_03$names1 <- unlist(map(str_split(list_rownames, "_"),5))
count_data_sequins_norm_03$names2 <- unlist(map(str_split(list_rownames, "_"),4))
count_data_sequins_norm_03$Species <- ifelse(count_data_sequins_norm_03$names1 == "", count_data_sequins_norm_03$names2, count_data_sequins_norm_03$names1)
count_data_sequins_norm_03$names1 <- NULL
count_data_sequins_norm_03$names2 <- NULL
count_data_sequins_norm_03$Organism <- paste(count_data_sequins_norm_03$Genus, count_data_sequins_norm_03$Species)
count_data_sequins_norm_03$Species <- NULL
count_data_sequins_norm_03$Genus <- NULL
rownames(count_data_sequins_norm_03) <- count_data_sequins_norm_03$Organism
count_data_sequins_norm_03$Organism <- NULL
count_data_sequins_norm_04 <- round(count_data_sequins_norm_03,5)

# export normalised absolute abundance counts
write.csv(count_data_sequins_norm_04, file = "count_data_sequins.csv", sep = ";")

# transpose final data table
count_data_sequins_norm_05 <- data.frame(t(count_data_sequins_norm_04))

####################################################################################
# import meta data table
meta_data_RandomForrest <- read_excel("meta_data_04_RandomForrest.xlsx", 
                                      col_types = c("text", "text", "text", "text", "text", "text", "numeric", 
                                                    "text", "numeric", "text", "numeric", "numeric", "numeric", 
                                                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
meta_data_RandomForrest <- data.frame(meta_data_RandomForrest)
meta_data_RandomForrest_sort <- meta_data_RandomForrest[order(meta_data_RandomForrest$sample_name),]

####################################################################################
# prepare random forest model
count_data_sequins_sort <- count_data_sequins_norm_05[order(rownames(count_data_sequins_norm_05)),]
final_df_rf <- data.frame(cbind(meta_data_RandomForrest_sort, count_data_sequins_sort))
final_df_rf$sample_name <- NULL
final_df_rf$patient_id <- NULL
final_df_rf$patient_number <- NULL
final_df_rf$sample_type <- factor(final_df_rf$sample_type, levels=c("throat swab", "sputum"), labels = c(0,1))
final_df_rf$sex <- factor(final_df_rf$sex, levels=c("f", "m"), labels = c(0,1))
final_df_rf$mutation_type <- factor(final_df_rf$mutation_type, levels=c("F/MF", "F/F"), labels = c(0,1))
final_df_rf$psa_chronic <- factor(final_df_rf$psa_chronic, levels=c("0", "1", "2", "3"), labels = c(0,1,2,3))
final_df_rf$pre_kaftrio_modulator <- factor(final_df_rf$pre_kaftrio_modulator, levels=c("naiv", "Orkambi", "Symkevi"), labels = c(0,1,2))
final_df_rf$time_point <- factor(final_df_rf$time_point, levels=c("V1", "V2", "V3"), labels = c(0,1,2))
final_df_rf$SA_count <- NULL
final_df_rf$PA_count <- NULL

####################################################################################
# run random forest
set.seed(1)

imp_final = NULL
error_rate_final = NULL

final_rf <- randomForest(time_point ~ ., data=final_df_rf, na.action = na.roughfix, ntree=80, mtry=12, importance=TRUE)
as.data.frame(final_rf$err.rate)

for (i in random_seeds){
  set.seed(i)
  final_rf <- randomForest(time_point ~ ., data=final_df_rf, na.action = na.roughfix, ntree=80, mtry=12, importance=TRUE)
  # store error rate locally
  error_rate_rf <- final_rf$err.rate
  # make data frame of error rate
  error_rate_rf <- as.data.frame(error_rate_rf)
  # add meta data
  error_rate_rf$Seed <- i
  # re-name columns
  colnames(error_rate_rf) <- c('OOB_all', 'Error_V1', 'Error_V2', 'Error_V3', 'Seed')
  error_rate_rf <- ddply(error_rate_rf, "Seed", numcolwise(mean))
  # transfer error rate data frame to global environment
  error_rate_final <- rbind(error_rate_final, error_rate_rf)
  
  # confirm random forest with boruta
  # run boruta
  final_df_rf_naFix <- na.roughfix(final_df_rf)
  boruta_rf <- Boruta(time_point~., data = final_df_rf_naFix, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_rf <- as.data.frame(boruta_rf$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_a3 <- importance(final_rf)
  imp_osps_a3 <- data.frame(predictors = rownames(imp_osps_a3), imp_osps_a3)
  
  # Make a data frame with predictor names and their importance
  imp_rf <- importance(final_rf)
  imp_rf <- data.frame(predictors = rownames(imp_rf), imp_rf)
  # add meta data
  imp_rf$Boruta_name <- rownames(boruta_rf)
  imp_rf$Boruta_predict <- boruta_rf$`boruta_rf$finalDecision`
  imp_rf_sub <- subset(imp_rf, MeanDecreaseAccuracy > 0.0)
  # add meta data
  imp_rf_sub$Seed <- i
  # store table globally
  imp_final <- rbind(imp_final, imp_rf_sub)
}

imp_final_conf <- imp_final[imp_final$Boruta_predict=="Confirmed",]
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "sweat_chloride", "Sweat chloride")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "sermet", "Sermet score")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "sample_type", "Sample type")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "lci", "Lung clearance index")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "total_bacterial_load", "Total bacterial load")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "sweat_rate", "Sweat rate")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Rothia.mucilaginosa", "Rothia mucilaginosa")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Staphylococcus.aureus", "Staphylococcus aureus")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Prevotella.melaninogenica", "Prevotella melaninogenica")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Rothia.aeria", "Rothia aeria")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "mef25", "Rothia aeria")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Rothia.dentocariosa", "Rothia dentocariosa")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "evenness", "Pielou eveness")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Prevotella.jejuni", "Prevotella jejuni")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Streptococcus.salivarius", "Streptococcus salivarius")
imp_final_conf$Boruta_name <- str_replace(imp_final_conf$Boruta_name, "Streptococcus.parasanguinis", "Streptococcus parasanguinis")

# plot and evaluate random forest results
rf_acc_plot <-
  ggplot(imp_final_conf, aes(y=reorder(Boruta_name,MeanDecreaseAccuracy), x=MeanDecreaseAccuracy)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(size=0.2) + theme_bw() +
  scale_y_discrete(label=c(expression(italic("Streptococcus parasanguinis")), 
                           expression(italic("Streptococcus salivarius")), 
                           expression(italic("Prevotella jejuni")), 
                           "Pielou's evenness index", 
                           expression(italic("Rothia dentocariosa")),
                           expression(italic("Rothia aeria")), 
                           expression(italic("Prevotella melaninogenica")), 
                           expression(italic("Staphylococcus aureus")), 
                           expression(italic("Rothia mucilaginosa")), 
                           "ß-adrenergic sweat rate", 
                           "Total bacterial load", 
                           "Lung clearance index",
                           "Sample type", 
                           "NPD Sermet score", 
                           "Sweat chloride concentration")) +
  xlim(0,10) +
  theme(panel.grid = element_blank()) + ylab("")

rf_gini_plot <-
  ggplot(imp_final_conf, aes(y=reorder(Boruta_name,MeanDecreaseGini), x=MeanDecreaseGini)) +
  geom_boxplot(outlier.size = 0) +
  geom_point(size=0.2) + theme_bw() +
  scale_y_discrete(label=c(expression(italic("Staphylococcus aureus")), 
                           expression(italic("Streptococcus parasanguinis")), 
                           expression(italic("Rothia dentocariosa")), 
                           expression(italic("Streptococcus salivarius")), 
                           expression(italic("Rothia aeria")), 
                           expression(italic("Prevotella jejuni")), 
                           "Pielou's evenness index", 
                           expression(italic("Prevotella melaninogenica")), 
                           "ß-adrenergic sweat rate",
                           expression(italic("Rothia mucilaginosa")), 
                           "Sample type",
                           "Lung clearance index",
                           "Total bacterial load", 
                           "NPD Sermet score", 
                           "Sweat chloride concentration")) +
  xlim(0,10) +
  theme(panel.grid = element_blank()) + ylab("")

error_rate_final_L <- gather(error_rate_final, key="OOB", value="value", -c("Seed"))
error_rate_final_L$OOB <- str_replace(error_rate_final_L$OOB, "Error_V1", "Baseline")
error_rate_final_L$OOB <- str_replace(error_rate_final_L$OOB, "Error_V2", "V2")
error_rate_final_L$OOB <- str_replace(error_rate_final_L$OOB, "Error_V3", "V3")
error_rate_final_L$OOB <- str_replace(error_rate_final_L$OOB, "OOB_all", "Total")
error_rate_final_L$OOB <- factor(error_rate_final_L$OOB, levels = c("Total", "Baseline", "V2", "V3"))
error_rate_med <- ddply(error_rate_final_L, "OOB", numcolwise(sd))

error_plot <-
  ggplot(error_rate_final_L) +
  geom_boxplot(aes(x=OOB, y=value), width=0.6, outlier.alpha = 0) +
  geom_jitter(aes(x=OOB, y=value), width=0.2, size=0.5) +
  ylim(0,1) + theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none") + 
  ylab("Out-of-bag error") + xlab("")

gini_acc_plot <- ggarrange(rf_acc_plot, rf_gini_plot, nrow=2, labels = c("A", "B"))
rf_final_plot <- ggarrange(gini_acc_plot, error_plot, nrow=1, widths = c(1,0.5), labels = c("A", "C"))

# export random forest figure
pdf("random_forest.pdf", width = 8, height = 6)
rf_final_plot
dev.off()

####################################################################################

# Alpha and beta diversity evaluation
# Prepare input table
count_data_div <- count_data_sequins_norm_05
count_data_div$bacterial_load <- rowSums(count_data_div)
count_data_div$patient
count_data_div$shannon_div <- vegan::diversity(count_data_div[,1:(ncol(count_data_div)-1)], index = "shannon")
count_data_div$simpson_div <- vegan::diversity(count_data_div[,1:(ncol(count_data_div)-2)], index = "simpson")
count_data_div$spec_richness <- vegan::specnumber(count_data_div[,1:(ncol(count_data_div)-3)])
count_data_div$sample_id <- rownames(count_data_div)
count_data_div$sample_type <- unlist(map(str_split(count_data_div$sample_id, "_"),1))
count_data_div$patient <- unlist(map(str_split(count_data_div$sample_id, "_"),3))
count_data_div$birth_year <- unlist(map(str_split(count_data_div$sample_id, "_"),4))
count_data_div$treatment <- unlist(map(str_split(count_data_div$sample_id, "_"),5))

# Extract important features
count_data_div <- select(count_data_div, c("patient","treatment", "sample_id", "shannon_div", "simpson_div", "bacterial_load", "spec_richness"))

# make column with evenness information
count_data_div$eveness <- count_data_div$shannon_div / log(count_data_div$spec_richness)

# re-order data table
count_data_div <- count_data_div[order(count_data_div$patient),]

# remove sputum samples
count_data_div <- count_data_div[!grepl("SP", count_data_div$sample_id),]

# split data frames based on sampling time point and keep only patients with longitudinal samples
count_data_div_v1 <- count_data_div[count_data_div$treatment=="V1",]
count_data_div_v2 <- count_data_div[count_data_div$treatment=="V2",]
count_data_div_v3 <- count_data_div[count_data_div$treatment=="V3",]

remove_1 <- setdiff(count_data_div_v1$patient, count_data_div_v2$patient)
count_data_div_v1  <- count_data_div_v1[!count_data_div_v1$patient %in% remove_1,]

remove_2 <- setdiff(count_data_div_v2$patient, count_data_div_v1$patient)
count_data_div_v2  <- count_data_div_v2[!count_data_div_v2$patient %in% remove_2,]

remove_3 <- setdiff(count_data_div_v3$patient, count_data_div_v1$patient)
count_data_div_v3  <- count_data_div_v3[!count_data_div_v3$patient %in% remove_2,]

df_v1_v2_div <- data.frame(Baseline=count_data_div_v1$shannon_div, V2=count_data_div_v2$shannon_div, V3=count_data_div_v3$shannon_div)
df_v1_v2_load <- data.frame(Baseline=count_data_div_v1$bacterial_load, V2=count_data_div_v2$bacterial_load, V3=count_data_div_v3$bacterial_load)
df_v1_v2_evenness <- data.frame(Baseline=count_data_div_v1$eveness, V2=count_data_div_v2$eveness, V3=count_data_div_v3$eveness)
df_v1_v2_simpson <- data.frame(Baseline=count_data_div_v1$simpson_div, V2=count_data_div_v2$simpson_div, V3=count_data_div_v3$simpson_div)

# Shannon diversity comparison (pair-wise): baseline to V2
shannon_div_V1_V2 <-
  ggpaired(df_v1_v2_div, "Baseline", "V2",
           palette = "jco", 
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=3) + 
  theme_bw() + xlab("Time point") + ylab("Shannon diversity index") +
  theme(panel.grid = element_blank()) + ylim(0,3.5) 

# Shannon diversity comparison (pair-wise): baseline to V3
shannon_div_V1_V3 <-
  ggpaired(df_v1_v2_div, "Baseline", "V3",
           palette = "jco", 
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=3) +  
  theme_bw() + xlab("Time point") + ylab(" ") +
  theme(panel.grid = element_blank()) + ylim(0,3.5)

# Merge both plots
shannon_div <- ggarrange(shannon_div_V1_V2, shannon_div_V1_V3)


# Simpson diversity comparison (pair-wise): baseline to V2
simpson_div_V1_V2 <-
  ggpaired(df_v1_v2_simpson, "Baseline", "V2",
           palette = "jco", 
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=3) + 
  theme_bw() + xlab("Time point") + ylab("Simpson diversity index") +
  theme(panel.grid = element_blank()) + ylim(0,1)

# Simpson diversity comparison (pair-wise): baseline to V3
simpson_div_V1_V3 <-
  ggpaired(df_v1_v2_simpson, "Baseline", "V3",
           palette = "jco", 
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=3) +  
  theme_bw() + xlab("Time point") + ylab(" ") +
  theme(panel.grid = element_blank()) + ylim(0,1)

# Merge both plots
simpson_div <- ggarrange(simpson_div_V1_V2, simpson_div_V1_V3)

# Bacterial load
df_v1_v2_load$Baseline_log <- log10(df_v1_v2_load$Baseline+1)
df_v1_v2_load$V2_log <- log10(df_v1_v2_load$V2+1)
df_v1_v2_load$V3_log <- log10(df_v1_v2_load$V3+1)

# Bacterial load comparison (pair-wise): baseline to V2
load_V1_V2 <-
  ggpaired(df_v1_v2_load, "Baseline_log", "V2_log",
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=3) +  
  theme_bw() + xlab(" ") + ylab("Bacterial load") +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(breaks = c(0,1,2,3,4), 
                     labels=c("0", expression("10"^1), expression("10"^2), expression("10"^3), expression("10"^4))) +
  scale_x_discrete("Time point",labels=c("Baseline", "V2"))

# Bacterial load comparison (pair-wise): baseline to V3
load_V1_V3 <-
  ggpaired(df_v1_v2_load, "Baseline_log", "V3_log",
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=4, colour="darkred") +  
  theme_bw() + xlab(" ") + ylab(" ") +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(breaks = c(0,1,2,3,4), 
                     labels=c("0", expression("10"^1), expression("10"^2), expression("10"^3), expression("10"^4))) +
  scale_x_discrete("Time point", labels=c("Baseline", "V3"))

# Merge both plots
bac_load <- ggarrange(load_V1_V2, load_V1_V3)
df_baseline_load <- data.frame(value=df_v1_v2_load$Baseline, group="Baseline")
# Calculate effect sizes and confidence intervals
df_v3_load <- data.frame(value=df_v1_v2_load$V2, group="V3")
load_merged <- data.frame(rbind(df_baseline_load, df_v3_load))
rcompanion::wilcoxonR(load_merged$value, g = load_merged$group, ci= TRUE) 

# Evenness comparison (pair-wise): baseline to V2
evenness_V1_V2 <-
  ggpaired(df_v1_v2_evenness, "Baseline", "V2",
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=3) + 
  theme_bw() + xlab(" ") + ylab("Pielou's evenness index") +
  theme(panel.grid = element_blank())  +
  scale_x_discrete("Time point", labels=c("Baseline", "V2")) + ylim(0,1)

# Evenness comparison (pair-wise): baseline to V3
evenness_V1_V3 <-
  ggpaired(df_v1_v2_evenness, "Baseline", "V3",
           line.color = "gray", line.size = 0.4,
           short.panel.labs = FALSE) +
  stat_compare_means(label = "p.signif", paired = TRUE, label.x.npc = "centre", size=4, colour="darkred") +  
  theme_bw() + xlab(" ") + ylab(" ") +
  theme(panel.grid = element_blank()) +
  scale_x_discrete("Time point", labels=c("Baseline", "V3")) + ylim(0,1)
# merge both plots
bac_evenness <- ggarrange(evenness_V1_V2, evenness_V1_V3)

# Calculate effect sizes and confidence intervals
df_baseline <- data.frame(value=df_v1_v2_evenness$Baseline, group="Baseline")
df_v3 <- data.frame(value=df_v1_v2_evenness$V3, group="V3")
evenness_merged <- data.frame(rbind(df_baseline, df_v3))
rcompanion::wilcoxonR(evenness_merged$value, g = evenness_merged$group, ci= TRUE) 

# merge plots for Figure 04
alpha_div_plots <- ggarrange(bac_load, bac_evenness, nrow=1, labels = c("A", "B"))
# merge plots for Supplementary Figure 02
alpha_div_supp <- ggarrange(shannon_div, simpson_div, nrow=1, labels=c("A", "B"))

####################################################################################
# Estimation of beta diversity
# Prepare input table and just keep cough swabs
count_data_sequins_norm_04_RA <- count_data_sequins_norm_04[,grepl("RA",colnames(count_data_sequins_norm_04))]
group_list <- unlist(map(str_split(colnames(count_data_sequins_norm_04_RA), "_"),5))

count_data_sequins_norm_04_RA_dist <- dist(t(log10(count_data_sequins_norm_04_RA+1)), method = "euclidean")
# run betadisper
betaplot <- betadisper(count_data_sequins_norm_04_RA_dist, group_list, type = "centroid", sqrt.dist = TRUE, bias.adjust=TRUE)

# store centroid in data frame
betaplot_df_centroid <- data.frame(betaplot$centroids)
betaplot_df_centroid$time_point <- rownames(betaplot_df_centroid)

# store positions of points in multivariate space in data frame
betaplot_df_points <- data.frame(betaplot$vectors)
betaplot_df_points$time_point <- unlist(map(str_split(rownames(betaplot_df_points), "_"),5))

# store distances in data frame
betaplot_df_dist <- data.frame(betaplot$distances)
betaplot_df_dist$time_point <- unlist(map(str_split(rownames(betaplot_df_dist), "_"),5))


time_point_comparison <- list(c("V1", "V2"),
                              c("V1", "V3"),
                              c("V2", "V3"))
# make distance boxplot
distance_centroid_plot <-
  ggplot(betaplot_df_dist, aes(x=time_point, y=betaplot.distances)) +
  geom_violin() +
  geom_jitter(width=0.02, size=1) +
  stat_compare_means() + 
  stat_compare_means(comparisons = time_point_comparison, label = "p.signif") +
  geom_pointrange(mapping = aes(x = time_point, y = betaplot.distances),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)}, # median and  25% quartile and 75% quartile
                  fun = median,
                  colour="red", size=0.5) +
  scale_x_discrete("Time point", labels=c("Baseline", "V2", "V3")) +
  scale_y_continuous("Distance to centroid", 
                     breaks = c(0, 1.5,3),
                     limits = c(0,4)) +
  theme_bw() + theme(panel.grid = element_blank())

# Calculate effect size of disperson
rcompanion::epsilonSquared(betaplot_df_dist$betaplot.distances, g=betaplot_df_dist$time_point, ci=TRUE)

# make beta diversity plot
betadisp_plot <-
  ggplot() +
  stat_ellipse(data=betaplot_df_points, aes(x=PCoA1, y=PCoA2, colour=time_point), type = "norm") +
  geom_point(data=betaplot_df_points, aes(x=PCoA1, y=PCoA2, colour=time_point, shape=time_point), size=2.5) + 
  geom_point(data=betaplot_df_centroid, aes(x=PCoA1, y=PCoA2, colour=time_point, shape=time_point), size=6) +
  scale_colour_manual(values = c("V1"="firebrick1", "V2"="dodgerblue3", "V3"="black"), labels = c("V1"="Baseline", "V2"="V2", "V3"="V3")) +
  scale_shape_manual(values = c("V1"=18, "V2"=15, "V3"=19), labels = c("V1"="Baseline", "V2"="V2", "V3"="V3")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.title = element_blank())

# combine beta diversity and distance plot
beta_div_test <- ggarrange(betadisp_plot, distance_centroid_plot, labels = c("D", "E"), widths = c(1,0.7))
# combine alpha and beta diversity plots
alpha_beta_plots <- ggarrange(alpha_div_plots, beta_div_test, nrow=2, heights = c(0.8,1))

# export supplementary Figure (Shannon and Simpson diversity indices)
pdf("Supplementary_Figure_S02.pdf", width = 8, height = 4)
alpha_div_supp
dev.off()

# export Figure for main test
pdf("Figure_04.pdf", width = 9, height = 8)
alpha_beta_plots
dev.off()


####################################################################################
# Generate species co-occurence network files
# select V1 time point
count_data_sequins_norm_V1 <- count_data_sequins_norm_04[,grepl("V1",colnames(count_data_sequins_norm_04))]
count_data_sequins_norm_V1_2 <- count_data_sequins_norm_V1[,grepl("RA",colnames(count_data_sequins_norm_V1))]
V1_mean_abundance <- rowMedians(as.matrix(count_data_sequins_norm_V1_2))

# perform Spearman's rank correlation analysis
spearman_norm_V1 <- rcorr(as.matrix(t(count_data_sequins_norm_V1_2)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_V1 <- spearman_norm_V1$P
# extract and store correlation coefficients of analysis
spearman_r_norm_V1 <- spearman_norm_V1$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_V1_edges <- reshape2::melt(spearman_p_norm_V1)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_V1_cor_edges <- reshape2::melt(spearman_r_norm_V1)
# merge tables
spearman_V1_edges$COR <- spearman_V1_cor_edges$value
# round p-values to five decimal places
spearman_V1_edges$value <- round(spearman_V1_edges$value, 5)
# store row names in own column
spearman_V1_edges$Label <- row.names(spearman_V1_cor_edges)
# rename columns
colnames(spearman_V1_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_V1_edges) <- NULL
# remove NAs
spearman_V1_edges <- spearman_V1_edges[complete.cases(spearman_V1_edges), ]
# extract significant correlations
spearman_V1_edges_sig <- spearman_V1_edges[spearman_V1_edges$pValue < 0.01,]
spearman_V1_edges_sig <- subset(spearman_V1_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_V1_edges_sig$Correlation <- ifelse(spearman_V1_edges_sig$Weight > 0, "pos", "neg")
spearman_V1_edges_sig$Type <- "Directed"
spearman_V1_edges_sig_final <- select(spearman_V1_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_V1_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_V1_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_V1_nodes_sig <- select(spearman_V1_edges_sig, c("Target"))
spearman_V1_nodes_sig$Target2 <- spearman_V1_nodes_sig$Target
spearman_V1_nodes_sig$Correlation <- spearman_V1_edges_sig$Correlation
spearman_V1_nodes_sig
# remove duplicate entries
spearman_V1_nodes_sig = spearman_V1_nodes_sig[!duplicated(spearman_V1_nodes_sig$Target),]
# make data frame
spearman_V1_nodes_sig <- data.frame(spearman_V1_nodes_sig)
# re-index rows
rownames(spearman_V1_nodes_sig) <- NULL
# rename columns
colnames(spearman_V1_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_V1_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_V1_nodes_sig$Label)," "), `[`, 1)
mean_list_V1 = c()
for (items in spearman_V1_nodes_sig$Id){
  mean_list_V1 = append(mean_list_V1, V1_mean_abundance[items])}
spearman_V1_nodes_sig$mean_abundance <- mean_list_V1
# export node list
write.table(spearman_V1_nodes_sig, file="nodes_V1.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_V1_edges_sig_final, file="edges_V1.csv", sep=";", col.names = TRUE, row.names = FALSE)


# select V2 time point
count_data_sequins_norm_V2 <- count_data_sequins_norm_04[,grepl("V2",colnames(count_data_sequins_norm_04))]
V2_mean_abundance <- rowMedians(as.matrix(count_data_sequins_norm_V2))
# perform Spearman's rank correlation analysis
spearman_norm_V2 <- rcorr(as.matrix(t(count_data_sequins_norm_V2)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_V2 <- spearman_norm_V2$P
# extract and store correlation coefficients of analysis
spearman_r_norm_V2 <- spearman_norm_V2$r
# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_V2_edges <- reshape2::melt(spearman_p_norm_V2)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_V2_cor_edges <- reshape2::melt(spearman_r_norm_V2)
# merge tables
spearman_V2_edges$COR <- spearman_V2_cor_edges$value
# round p-values to five decimal places
spearman_V2_edges$value <- round(spearman_V2_edges$value, 5)
# store row names in own column
spearman_V2_edges$Label <- row.names(spearman_V2_cor_edges)
# rename columns
colnames(spearman_V2_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_V2_edges) <- NULL
# remove NAs
spearman_V2_edges <- spearman_V2_edges[complete.cases(spearman_V2_edges), ]
# extract significant correlations
spearman_V2_edges_sig <- spearman_V2_edges[spearman_V2_edges$pValue < 0.01,]
spearman_V2_edges_sig <- subset(spearman_V2_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_V2_edges_sig$Correlation <- ifelse(spearman_V2_edges_sig$Weight > 0, "pos", "neg")
spearman_V2_edges_sig$Type <- "Directed"
spearman_V2_edges_sig_final <- select(spearman_V2_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_V2_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_V2_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_V2_nodes_sig <- select(spearman_V2_edges_sig, c("Target"))
spearman_V2_nodes_sig$Target2 <- spearman_V2_nodes_sig$Target
spearman_V2_nodes_sig$Correlation <- spearman_V2_edges_sig$Correlation
spearman_V2_nodes_sig
# remove duplicate entries
spearman_V2_nodes_sig = spearman_V2_nodes_sig[!duplicated(spearman_V2_nodes_sig$Target),]
# make data frame
spearman_V2_nodes_sig <- data.frame(spearman_V2_nodes_sig)
# re-index rows
rownames(spearman_V2_nodes_sig) <- NULL

# rename columns
colnames(spearman_V2_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_V2_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_V2_nodes_sig$Label)," "), `[`, 1)

mean_list_V2 = c()
for (items in spearman_V2_nodes_sig$Id){
  mean_list_V2 = append(mean_list_V2, V2_mean_abundance[items])}
spearman_V2_nodes_sig$mean_abundance <- mean_list_V2
# export node list
write.table(spearman_V2_nodes_sig, file="nodes_V2.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_V2_edges_sig_final, file="edges_V2.csv", sep=";", col.names = TRUE, row.names = FALSE)


# select V3 time point
count_data_sequins_norm_V3 <- count_data_sequins_norm_04[,grepl("V3",colnames(count_data_sequins_norm_04))]
V3_mean_abundance <- rowMedians(as.matrix(count_data_sequins_norm_V3))
# perform Spearman's rank correlation analysis
spearman_norm_V3 <- rcorr(as.matrix(t(count_data_sequins_norm_V3)), type = 'spearman')
# create node and edge lists
# extract and store p-values of correlation analysis
spearman_p_norm_V3 <- spearman_norm_V3$P
# extract and store correlation coefficients of analysis
spearman_r_norm_V3 <- spearman_norm_V3$r
# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
spearman_V3_edges <- reshape2::melt(spearman_p_norm_V3)
# take data in wide format and stack columns into a single column of data for correlation coefficients
spearman_V3_cor_edges <- reshape2::melt(spearman_r_norm_V3)
# merge tables
spearman_V3_edges$COR <- spearman_V3_cor_edges$value
# round p-values to five decimal places
spearman_V3_edges$value <- round(spearman_V3_edges$value, 5)
# store row names in own column
spearman_V3_edges$Label <- row.names(spearman_V3_cor_edges)
# rename columns
colnames(spearman_V3_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(spearman_V3_edges) <- NULL
# remove NAs
spearman_V3_edges <- spearman_V3_edges[complete.cases(spearman_V3_edges), ]
# extract significant correlations
spearman_V3_edges_sig <- spearman_V3_edges[spearman_V3_edges$pValue < 0.01,]
spearman_V3_edges_sig <- subset(spearman_V3_edges_sig, Weight <= -0.2 | Weight >= 0.2)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
spearman_V3_edges_sig$Correlation <- ifelse(spearman_V3_edges_sig$Weight > 0, "pos", "neg")
spearman_V3_edges_sig$Type <- "Directed"
spearman_V3_edges_sig_final <- select(spearman_V3_edges_sig, c("Source", "Target", "Type", "Correlation"))
spearman_V3_edges_sig_final$Genus <- sapply(strsplit(as.character(spearman_V3_edges_sig_final$Target)," "), `[`, 1)
#select both columns to create node list
spearman_V3_nodes_sig <- select(spearman_V3_edges_sig, c("Target"))
spearman_V3_nodes_sig$Target2 <- spearman_V3_nodes_sig$Target
spearman_V3_nodes_sig$Correlation <- spearman_V3_edges_sig$Correlation
spearman_V3_nodes_sig
# remove duplicate entries
spearman_V3_nodes_sig = spearman_V3_nodes_sig[!duplicated(spearman_V3_nodes_sig$Target),]
# make data frame
spearman_V3_nodes_sig <- data.frame(spearman_V3_nodes_sig)
# re-index rows
rownames(spearman_V3_nodes_sig) <- NULL
# rename columns
colnames(spearman_V3_nodes_sig) <- c("Id", "Label", "Correlation")
# add genus information
spearman_V3_nodes_sig$Genus <- sapply(strsplit(as.character(spearman_V3_nodes_sig$Label)," "), `[`, 1)
mean_list_V3 = c()
for (items in spearman_V3_nodes_sig$Id){
  mean_list_V3 = append(mean_list_V3, V3_mean_abundance[items])}
spearman_V3_nodes_sig$mean_abundance <- mean_list_V3
# export node list
write.table(spearman_V3_nodes_sig, file="nodes_V3.csv", sep=";", col.names = TRUE, row.names = FALSE)
# export edge list
write.table(spearman_V3_edges_sig_final, file="edges_V3.csv", sep=";", col.names = TRUE, row.names = FALSE)


####################################################################################
# Import gephi output
V1_gephi_output <- read_delim("V1_gephi_output.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
V1_gephi_output <- data.frame(V1_gephi_output)
V1_gephi_output$V <- "V1"
V1_gephi_output_largestComponent <- V1_gephi_output[V1_gephi_output$modularity_class==3,]

V2_gephi_output <- read_delim("V2_gephi_output.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
V2_gephi_output <- data.frame(V2_gephi_output)
V2_gephi_output$V <- "V2"
V2_gephi_output$Cluster <- NULL
V2_gephi_output_largestComponent <- V2_gephi_output[V2_gephi_output$modularity_class==0,]

V3_gephi_output <- read_delim("V3_gephi_output.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
V3_gephi_output <- data.frame(V3_gephi_output)
V3_gephi_output$V <- "V3"
V3_gephi_output$Cluster <- NULL
V3_gephi_output_largestComponent <- V3_gephi_output[V3_gephi_output$modularity_class==2,]

gephi_output <- data.frame(rbind(V1_gephi_output, V2_gephi_output, V3_gephi_output))

V1_largest_sum <- sum(V1_gephi_output_largestComponent$mean_abundance)
V1_gephi_output_largestComponent$mean_abundance_percentage <- (V1_gephi_output_largestComponent$mean_abundance / V1_largest_sum) * 100
V1_sum_authority <- sum(V1_gephi_output_largestComponent$Authority)
V1_gephi_output_largestComponent$Authority_rel <- V1_gephi_output_largestComponent$Authority / V1_sum_authority
V1_gephi_output_largestComponent_spec <- V1_gephi_output_largestComponent
V1_gephi_output_largestComponent <- ddply(V1_gephi_output_largestComponent, "genus", numcolwise(sum))
V1_gephi_output_largestComponent$Treatment <- "V1"
V1_gephi_output_largestComponent_spec$Treatment <- "V1"

V2_largest_sum <- sum(V2_gephi_output_largestComponent$mean_abundance)
V2_gephi_output_largestComponent$mean_abundance_percentage <- (V2_gephi_output_largestComponent$mean_abundance / V2_largest_sum) * 100
V2_sum_authority <- sum(V2_gephi_output_largestComponent$Authority)
V2_gephi_output_largestComponent$Authority_rel <- V2_gephi_output_largestComponent$Authority / V2_sum_authority
V2_gephi_output_largestComponent_spec <- V2_gephi_output_largestComponent
V2_gephi_output_largestComponent <- ddply(V2_gephi_output_largestComponent, "genus", numcolwise(sum))
V2_gephi_output_largestComponent$Treatment <- "V2"
V2_gephi_output_largestComponent_spec$Treatment <- "V2"

V3_largest_sum <- sum(V3_gephi_output_largestComponent$mean_abundance)
V3_gephi_output_largestComponent$mean_abundance_percentage <- (V3_gephi_output_largestComponent$mean_abundance / V3_largest_sum) * 100
V3_sum_authority <- sum(V3_gephi_output_largestComponent$Authority)
V3_gephi_output_largestComponent$Authority_rel <- V3_gephi_output_largestComponent$Authority / V3_sum_authority
V3_gephi_output_largestComponent_spec <- V3_gephi_output_largestComponent
V3_gephi_output_largestComponent <- ddply(V3_gephi_output_largestComponent, "genus", numcolwise(sum))
V3_gephi_output_largestComponent$Treatment <- "V3"
V3_gephi_output_largestComponent_spec$Treatment <- "V3"

# merge largest components from all time points
largest_component <- data.frame(rbind(V1_gephi_output_largestComponent, V2_gephi_output_largestComponent, V3_gephi_output_largestComponent))
largest_component$Component <- "1"

all_components <- data.frame(rbind(largest_component))
# clean data frame and remove low-abundance taxa in first network component
all_components$genus <- str_replace_all(all_components$genus, "Agarilytica", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Alkalitalea", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Bernardetia", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Butyrivibrio", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Chitinophaga", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Draconibacterium", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Dyadobacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Emticicia", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Filimonas", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Flammeovirgaceae", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Flavisolibacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Haliscomenobacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Hymenobacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Labilithrix", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Marinifilaceae", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Minicystis", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Moorea", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Runella", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Sebaldella", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Sphingobacteriaceae", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Spirosoma", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Selenomonas", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Porphyromonadaceae", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Pontibacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Pelosinus", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Pedobacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Paludibacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Blautia", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Nostocales", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Roseburia", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Paenibacillus", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Treponema", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Alistipes", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Anaerostipes", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Lachnoclostridium", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Niastella", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Megasphaera", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Abiotrophia", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Flavobacterium", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Megasphaera", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Odoribacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Parabacteroides", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Chryseobacterium", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Pseudopropionibacterium", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Barnesiella", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Pseudomonas", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Staphylococcus", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Lautropia", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Clostridium", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Campylobacter", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Streptobacillus", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Filifactor", "Other")
all_components$genus <- str_replace_all(all_components$genus, "Aggregatibacter", "Other")
all_components_V1 <- all_components[all_components$Component=="1",]
all_components_V2 <- all_components[all_components$Treatment=="V3",]

abundance_plot <-
  ggplot(all_components_V1, aes(x=Treatment, y=mean_abundance_percentage, fill = genus)) +
  geom_col(color = "black", position = position_stack(reverse = FALSE), width=0.8) + 
  theme_bw() + 
  scale_fill_manual(values=c("Actinomyces"="cyan2",
                             "Atopobium"="grey20",             
                             "Bacteroides"="aquamarine4",
                             "Bifidobacterium"="#9e0142",
                             "Capnocytophaga"="#d53e4f",
                             "Eubacterium"="green",
                             "Fusobacterium"="#f46d43",
                             "Gemella"="peru",
                             "Haemophilus"="#fdae61",
                             "Leptotrichia"="wheat1",
                             "Neisseria"="#fee08b",
                             "Parvimonas"="thistle3",
                             "Porphyromonas"="#e6f598",
                             "Prevotella"="#abdda4",
                             "Rothia"="#66c2a5",
                             "Schaalia"="yellow",
                             "Streptococcus"="#3288bd",
                             "Tannerella"="mediumorchid1",
                             "Veillonella"="#5e4fa2",
                             "Other"="beige")) + xlab("Time point") + ylab("Median abundance (in %)") +
  scale_x_discrete(labels=c("Baseline", "V2", "V3")) +
  guides(fill=guide_legend(ncol=1)) +
  theme(panel.grid = element_blank(), legend.position = "right", legend.title = element_blank(), 
        strip.background = element_rect(fill = "white")) 

# look into largest component on species level
largest_component_spec <- data.frame(rbind(V1_gephi_output_largestComponent_spec, 
                                           V2_gephi_output_largestComponent_spec, 
                                           V3_gephi_output_largestComponent_spec))
largest_component_spec$Component <- "1"

compare_treatment <- list(c("V1", "V2"),
                          c("V2", "V3"),
                          c("V1", "V3"))

largest_component_spec_div <- select(largest_component_spec, c("Treatment", "Label", "mean_abundance"))
largest_component_spec_div_spread <- spread(largest_component_spec_div, key="Label", value="mean_abundance")
largest_component_spec_div_spread[is.na(largest_component_spec_div_spread)] <- 0
rownames(largest_component_spec_div_spread) <- largest_component_spec_div_spread$Treatment
largest_component_spec_div_spread$Treatment <- NULL
largest_component_spec_div_spread_t <- data.frame(t(largest_component_spec_div_spread))

# calculate Shannon diversity of the largest network component
shannon_div <- vegan::diversity(largest_component_spec_div_spread, index = "shannon")
div <- round(shannon_div,2)
div_df <- data.frame(div)
div_df$Treatment <- rownames(div_df)

# generate plot for node authority
spec_authority <-
  ggplot(largest_component_spec, aes(x=Treatment, y=Authority)) +
  geom_violin(width=1.3) + 
  geom_jitter(width=0.05) +
  stat_summary(fun="median", geom="point", colour="red", shape=15, size=2) +
  scale_y_continuous("Node authority", breaks = c(0.00, 0.08, 0.18), limits = c(0,0.23)) +
  #stat_compare_means() +
  stat_compare_means(comparisons = compare_treatment, method = "wilcox.test", 
                     label = "p.signif", label.y = c(0.15, 0.17, 0.19)) +
  theme_bw() +
  scale_x_discrete("Time point",labels=c("Baseline", "V2", "V3")) +
  theme(panel.grid = element_blank()) 

# generate plot for betweeness centrality
betweeness_cent <-
  ggplot(largest_component_spec, aes(x=Treatment, y=betweenesscentrality)) +
  geom_violin(width=1) + 
  geom_jitter(width=0.05) +
  stat_summary(fun="median", geom="point", colour="red", shape=15, size=2) +
  scale_y_continuous("Betweeness centrality", breaks = c(0.00, 0.04, 0.08), limits = c(0,0.1)) +
  # stat_compare_means() + 
  stat_compare_means(comparisons = compare_treatment, method = "wilcox.test", 
                     label = "p.signif", label.y = c(0.04, 0.07, 0.08)) +
  theme_bw() +
  scale_x_discrete("Time point",labels=c("Baseline", "V2", "V3")) +
  theme(panel.grid = element_blank())

# obtain effect sizes and confidence intervals
rcompanion::epsilonSquared(largest_component_spec$betweenesscentrality, g=largest_component_spec$Treatment, ci=TRUE) 
rcompanion::epsilonSquared(largest_component_spec$Authority, g=largest_component_spec$Treatment, ci=TRUE) 

# generate plot for Shannon diversity visualisation
shannon_div <-
  ggplot(div_df, aes(x=Treatment, y=div)) +
  geom_col() +
  ylab("Shannon diversity index") +
  theme_bw() +
  scale_x_discrete("Time point",labels=c("Baseline", "V2", "V3")) +
  theme(panel.grid = element_blank()) + coord_flip()

# generate an empty plot for network insertions
empty_plot <-
  ggplot(div_df, aes(x=Treatment, y=div)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "white"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) 

# merge all plots
empty_plots <- ggarrange(empty_plot, empty_plot, empty_plot, widths = c(1,0.8,0.8), labels = c("A", "B", "C"), nrow=1)
box_plots_two <- ggarrange(spec_authority, betweeness_cent, nrow=1, labels=c("F", "G"))
shannon_box_plots <- ggarrange(shannon_div, box_plots_two, nrow=2, heights = c(0.4,1), labels = c("E", "F"))
all_merged <- ggarrange(abundance_plot, shannon_box_plots, nrow=1, widths = c(0.4, 0.5), labels=c("D", "E"))
all_merged_2 <- ggarrange(empty_plots, all_merged, nrow=2)

pdf("Figure_05.pdf", width = 12, height = 12)
all_merged_2
dev.off()
