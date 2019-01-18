#### Script to look at both specices data and metadata break downs using RF


#load in metadata
metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17.csv",
                       sep="\t", header=T, row.names=1, comment.char="")
rownames(metadata) <- gsub("-",".",rownames(metadata))

metadata_fix <- metadata[which(metadata$Within_183_days=="Y"),]

#load in speices table
species <- read.table("~/projects/MALL-MGS-ALL/species_metaphlan2.tsv",
                      sep="\t", row.names=1, header=T, comment.char = "", quote="")

species_flip <- data.frame(t(species))

speices_fix <- species_flip[rownames(metadata_fix),]

rownames(speices_fix) == rownames(metadata_fix)

#add in metadata that should be included in ML model
speices_fix$Total_ABX <- metadata_fix$ABX_status
speices_fix$Infection <- metadata_fix$Infection_in_6
speices_fix$Vancomycin_exposure <- metadata_fix$Vanco_exposure
speices_fix$Fungal_expo <- metadata_fix$anti_fungal_exposure
speices_fix$Piptazo_expo <- metadata_fix$PIP_TAZ_exposure
speices_fix$Other_ABX <- metadata_fix$OTHER_ABX
speices_fix$Age <- metadata_fix$Age_at_diagnosis
speices_fix$Days_Since_Start_of_Treatment <- metadata_fix$Days_since_therapy_start
speices_fix$Treatment_Type <- metadata_fix$Treatment_type

#alright table is good to go lets go and go some RF

library(randomForest)
#set seed for reproducibility
set.seed(1995)
colnames(speices_fix)[299]

Infection_model <- randomForest(x=speices_fix[,-299], y=speices_fix[,299], importance = T, proximity = T, ntree=10001)
Infection_model
#84% accuracy lets see what the top hits were 

important_features <- data.frame(Infection_model$importance)



important_features <- important_features[order(important_features$MeanDecreaseAccuracy, decreasing = T),]
feature_names <- rownames(important_features)
feature_names
feature_names <- gsub(".*s__","",feature_names)
feature_names
feature_names <- gsub("_", " ", feature_names)


important_features$feature_name <- feature_names

library(forcats)
#set levels in order 
important_features$feature_name <- fct_inorder(important_features$feature_name) %>% fct_rev()

library(ggplot2)
library(cowplot)
#plot by rank of importance
important_features_plot <- ggplot(important_features[c(1:20),], aes(x=feature_name, y=MeanDecreaseAccuracy)) +
  geom_point() + 
  ylab("Meann Decrease in Accuracy") +
  xlab("Prediction Feature") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip()
important_features_plot






###okay lets do the same but keep no infection and vanco

#only keep NIC samples and IC samples that had Vanco Exposure
metadata_vanco_full <- metadata_fix[which(metadata_fix$Vanco_exposure=="Y" & metadata_fix$Infection_in_6=="Y" | metadata_fix$Infection_in_6=="No"),]
colnames(speices_fix)[300]
#remove Vancomycin Exposure from the model
speices_vanco_full <- speices_fix[rownames(metadata_vanco_full), -300]

vanco_expo <- randomForest(x=speices_vanco_full[,-299], y=speices_vanco_full[,299], importance=T, proximity = T, ntree=10001)
vanco_expo
View(vanco_expo$importance)

important_features_vanco <- data.frame(vanco_expo$importance)
important_features_vanco <- important_features_vanco[order(important_features_vanco$MeanDecreaseAccuracy, decreasing = T),]
feature_names_vanco <- rownames(important_features_vanco)
feature_names_vanco
feature_names_vanco <- gsub(".*s__","",feature_names_vanco)
feature_names_vanco
feature_names_vanco <- gsub("_", " ", feature_names_vanco)

important_features_vanco$feature_name <- feature_names_vanco
important_features_vanco$feature_name <- fct_inorder(important_features_vanco$feature_name) %>% fct_rev()

important_features_plot_vanco <- ggplot(important_features_vanco[c(1:20),], aes(x=feature_name, y=MeanDecreaseAccuracy)) +
  geom_point() + 
  ylab("Meann Decrease in Accuracy") +
  xlab("Prediction Feature") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip()
important_features_plot_vanco

#plot of features important in classifying samples based on NIC and IC with vanco exposure
#intresting that top feature is reduction of a gram neg!



#get NIC samples and IC samples with no Vanco Exposure
metadata_no_vanco_full <- metadata_fix[which(metadata_fix$Infection_in_6=="No" | (metadata_fix$Infection_in_6=="Y" & metadata_fix$Vanco_exposure=="N")),]
speices_no_vanco_full <- speices_fix[rownames(metadata_no_vanco_full), -300]

vanco_no_expo <- randomForest(x=speices_no_vanco_full[,-299], y=speices_no_vanco_full[,299], importance=T, proximity = T, ntree=10001)
vanco_no_expo
View(vanco_no_expo$importance)

important_features_no_vanco <- data.frame(vanco_no_expo$importance)
important_features_no_vanco <- important_features_no_vanco[order(important_features_no_vanco$MeanDecreaseAccuracy, decreasing = T),]
feature_names_no_vanco <- rownames(important_features_no_vanco)
feature_names_no_vanco
feature_names_no_vanco <- gsub(".*s__","",feature_names_no_vanco)
feature_names_no_vanco
feature_names_no_vanco <- gsub("_", " ", feature_names_no_vanco)

important_features_no_vanco$feature_name <- feature_names_no_vanco
important_features_no_vanco$feature_name <- fct_inorder(important_features_no_vanco$feature_name) %>% fct_rev()

important_features_plot_no_vanco <- ggplot(important_features_no_vanco[c(1:20),], aes(x=feature_name, y=MeanDecreaseAccuracy)) +
  geom_point() + 
  ylab("Meann Decrease in Accuracy") +
  xlab("Prediction Feature") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip()
important_features_plot_no_vanco

#make supplmental figure 7
sup_fig_7 <- plot_grid(important_features_plot_vanco, important_features_plot_no_vanco, labels="AUTO")
sup_fig_7

###they are the same indicating that vanco exposure doesn't matter?

#Make a Roc Curve (not in the manuscript)
library(ROCR)
predictions <- as.vector(Infection_model$votes[,2])
pred <- prediction(predictions, speices_fix$Infection)
perf_AUC <- performance(pred, "auc")
AUC=perf_AUC@y.values[[1]]
perf_ROC <- performance(pred, "tpr", "fpr")
plot(perf_ROC, main="ROC plot")
