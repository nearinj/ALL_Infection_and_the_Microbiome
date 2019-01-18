#beta div analysis

library(vegan)

#load in metadata
metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17.csv", sep="\t", 
                       header=T, row.names=1, comment.char="")
#reformat row naming

rownames(metadata) <- gsub("-", ".", rownames(metadata))

#load in weighted unifrac distance matrix
weighted_unifrac <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/diversity/weighted_unifrac/distance-matrix.tsv", 
                               sep="\t", header=T, row.names = 1)
rownames(weighted_unifrac) <- gsub("-", ".", rownames(weighted_unifrac))

#make row names match up with metadata

metadata_fix <- metadata[rownames(weighted_unifrac),]

#okay awesome! lets do some adnois tests
rownames(weighted_unifrac) == rownames(metadata_fix)

#### Do univariate PERMAONA analysis

mod1 <- adonis2(weighted_unifrac ~ metadata_fix$Infection_in_6)
mod1
#significant.

mod2 <- adonis2(weighted_unifrac ~ metadata_fix$Sex)
mod2
#significant.

mod3 <- adonis2(weighted_unifrac ~ metadata_fix$Age_at_diagnosis)
mod3


mod4 <- adonis2(weighted_unifrac ~ metadata_fix$Treatment_type)
mod4

mod5 <- adonis2(weighted_unifrac ~ metadata_fix$Days_since_therapy_start)
mod5



#test for exposure to antibiotics
mod6 <- adonis2(weighted_unifrac ~ metadata_fix$Vanco_exposure)
mod6

mod7 <- adonis2(weighted_unifrac ~ metadata_fix$anti_fungal_exposure)
mod7

mod8 <- adonis2(weighted_unifrac ~ metadata_fix$PIP_TAZ_exposure)
mod8

mod9 <- adonis2(weighted_unifrac ~ metadata_fix$OTHER_ABX)
mod9

#found no relationship between beta diversity and age in months... HOWEVER not that this is about more complicated and perhaps
###this should be done in a better way...
mod10 <- adonis2(weighted_unifrac ~ metadata_fix$Age_in_months)
mod10

#overall antibitoic exposure
mod10 <- adonis2(weighted_unifrac ~ metadata_fix$ABX_status)
mod10
#okay so significant interactions are as follows:
#Sex
#Days_since_therapy_start
#vanco exposure
#anti_fungal exposure
#Infection in 6 months


#multivariate PERMANOVA
sig_mod <- adonis2(weighted_unifrac ~ metadata_fix$Infection_in_6 + 
                     metadata_fix$Sex + 
                     metadata_fix$Days_since_therapy_start + 
                     metadata_fix$Vanco_exposure +
                     metadata_fix$anti_fungal_exposure, by="margin")
sig_mod
#okay lets plot this bad boy

#build PCoA
#use APE to make PCoA
library(ape)
library(ggplot2)
attach(metadata_fix)
all_pcoa <- pcoa(weighted_unifrac)
all_cords <- as.data.frame(all_pcoa$vectors)
metadata_fix$PC1 <- all_cords$Axis.1
metadata_fix$PC2 <- all_cords$Axis.2
metadata_fix <- metadata_fix[order(as.Date(metadata_fix$Date_of_sample)),]
metadata_fix$Patient <- gsub("P","",metadata_fix$Patient)
metadata_fix$Infection_in_6 <- ifelse(metadata_fix$Infection_in_6=="Y", "Yes", "No")
all_pcoa

#Fig 2A
infection_pcoa <- ggplot(metadata_fix, aes(PC1, PC2, color=metadata_fix$Infection_Col)) + geom_point(size=4) +
  ylab("PC2 (11.7%)") + xlab("PC1 (56.45%)") + guides(color=guide_legend(title="Infection Event")) + 
  theme_minimal() + scale_color_identity(guide = T, labels=c("No","Yes"))
infection_pcoa

#Fig 2B
all_pcoa_plot <- ggplot(metadata_fix, aes(PC1, PC2, color=metadata_fix$State_Cols))+ geom_point(size=5) + 
  guides(size=F, color=guide_legend(title="Infection Event")) +
  geom_path(aes(group=Patient), arrow = arrow(type="closed", angle=30, length=unit(0.4,"cm")), size=.5, color=metadata_fix$P_Color) + 
  geom_text(aes(label=Patient), color="black", size=1.75) + xlab("PC1 (56.45%)") + ylab("PC2 (11.07%)") +
  theme_minimal() + scale_color_identity(guide=T, labels=c("Never","Post","Pre")) 

all_pcoa_plot

#Figure2 
library(cowplot)
Figure2 <- plot_grid(infection_pcoa, all_pcoa_plot, labels="AUTO", rel_widths = 1.5)
Figure2

#sup fig with PC1 on Y and Days on X

beta_time_plot <- ggplot(metadata_fix, aes(Days_since_therapy_start, PC1, color=metadata_fix$State_Cols)) + geom_point(size=5) +
  guides(size=F, color=guide_legend(title = "Infection Event")) + scale_color_identity(guide=T, labels=c("Never","Post", "Pre")) +
  ylab("PC1 (56.45%)") + xlab("Days Since Therapy Start")
beta_time_plot

cor.test(metadata_fix$Days_since_therapy_start, metadata_fix$PC1, method='spearman')

beta_time_plot_infection <- ggplot(metadata_fix, aes(Days_since_therapy_start, PC1, color=metadata_fix$Infection_Col)) + geom_point(size=5) +
  guides(size=F, color=guide_legend(title="Infectious Complication")) + scale_color_identity(guide=T, labels=c("No", "Yes")) +
  ylab("PC1 (56.45%)") + xlab("Days Since Therapy Start")
beta_time_plot_infection

beta_age_plot <- ggplot(metadata_fix, aes(Age_at_diagnosis, PC1, color=metadata_fix$State_Cols)) + geom_point(size=5) +
  guides(size=F, color=guide_legend(title= "Infectious Event")) + scale_color_identity(guide=T, labels=c("Never", "Post", "Pre")) +
  ylab("PC1 (56.45%)") + xlab("Age at Diagnosis")
beta_age_plot

beta_age_plot_infection <- ggplot(metadata_fix, aes(Age_at_diagnosis, PC1, color=metadata_fix$Infection_Col)) + geom_point(size=5) +
  guides(size=F, color=guide_legend(title= "Infectious Complication")) + scale_color_identity(guide=T, labels=c("No", "Yes")) +
  ylab("PC1 (56.45%)") + xlab("Age at Diagnosis")
beta_age_plot_infection

sup_fig_age_beta <- plot_grid(beta_age_plot_infection, beta_time_plot_infection, labels="AUTO")
sup_fig_age_beta
#this plot looks good and can also explain why infection sits in 
#include this plot

wilcox.test(metadata_fix$PC1 ~ metadata_fix$Infection_in_6)


sup_plot_review <- plot_grid(beta_time_plot_infection, beta_time_plot, labels = "AUTO")
sup_plot_review
#do correlation of change in beta diversity for NIC and IC patients seperately.... (not sample replicaton per patient)
#can you correlate PC1 with time??? 


#stat analysis
hist(metadata_fix$Days_since_therapy_start)
#non-normal distrubition so we will use spearman correlation.

#no correlation between change in PC1 and and time...
cor.test(metadata_fix$PC1, metadata_fix$Days_since_therapy_start, method="spearman", exact=F)

#this however doesn't take into account the total variability of the samples over time...
#i.e samples may be moving up and then down 

#make sup fig 4
#vanco exposure pcoa
metadata_fix$Vanco_exposure <- ifelse(metadata_fix$Vanco_exposure=="Y", "Yes","No")
vanco_expo_plot <- ggplot(metadata_fix, aes(PC1, PC2, color=Vanco_exposure)) + geom_point(size=5) + theme_minimal() +
  guides(color=guide_legend(title = "Vancomycin Exposure")) + xlab("PC1 (56.45%)") + ylab("PC2 (11.07%)")
vanco_expo_plot

#antifungal exposure
metadata_fix$anti_fungal_exposure <- ifelse(metadata_fix$anti_fungal_exposure=="Y", "Yes", "No")
antifungal_expo_plot <- ggplot(metadata_fix, aes(PC1, PC2, color=anti_fungal_exposure)) + geom_point(size=5) + theme_minimal() +
  guides(color=guide_legend(title="Anti-Fungal Exposure")) + xlab("PC1 (56.45%)") + ylab("PC2 (11.07%)")
antifungal_expo_plot

Sup_Fig_4 <- plot_grid(vanco_expo_plot, antifungal_expo_plot, labels="AUTO")
Sup_Fig_4





#Lets do pre post never analysis

#pre v Post
metadata_pre_post <- droplevels(metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Pre" | metadata_fix$Pre_Post_never_infect=="Post"),])
weighted_unifrac_pre_post <- weighted_unifrac[rownames(metadata_pre_post), rownames(metadata_pre_post)]

adonis2(weighted_unifrac_pre_post ~ metadata_pre_post$Pre_Post_never_infect)
#no significant differents between pre and post based on weighted unifrac

#post v never
metadata_post_never <- droplevels(metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Post" | metadata_fix$Pre_Post_never_infect=="Never"),])
weighted_unifrac_post_never <- weighted_unifrac[rownames(metadata_post_never), rownames(metadata_post_never)]

adonis2(weighted_unifrac_post_never ~ metadata_post_never$Pre_Post_never_infect)
#significant differences between post and never based on weighted unifrac

#pre v never
metadata_pre_never <- droplevels(metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Pre" | metadata_fix$Pre_Post_never_infect=="Never"),])
weighted_unifrac_pre_never <- weighted_unifrac[rownames(metadata_pre_never), rownames(metadata_pre_never)]

adonis2(weighted_unifrac_pre_never ~ metadata_pre_never$Pre_Post_never_infect)
#done


#### Lets look at baseline samples for differences in sex. 

#load in baseline distance metrics
base_line_wuni <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/BASE_LINE/core-metric-phylo/weighted_unifrac/distance-matrix.tsv",
                             sep="\t", header=T, row.names = 1)
rownames(base_line_wuni) <- gsub("-",".",rownames(base_line_wuni))

metadata_base <- metadata[rownames(base_line_wuni),]

sex_diff <- adonis2(base_line_wuni ~ Sex, data=metadata_base, permutations=99999)
sex_diff
#no differences at baseline

Infec_diff <- adonis2(base_line_wuni ~ Infection_in_6, data=metadata_base, permuatations=99999)
Infec_diff
#no difference at baseline
rownames(base_line_wuni) == rownames(metadata_base)

#not super convincing as there was only 11 baseline samples for 16S sequencing.