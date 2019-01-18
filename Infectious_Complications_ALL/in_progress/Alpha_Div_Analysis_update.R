### Clean up all the previous scripts so we can upload them to github!


#load in patient metadata

metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17.csv", sep="\t",
                       row.names=1, header=T, comment.char="")

#load in faiths phylogenetic diversity scores
faiths_pd <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/diversity/alpha_data/faith_pd.tsv", 
                        sep="\t", row.names=1, header=T)
#load in shannon phylogenetic diversity scores
shannon <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/diversity/alpha_data/shannon_alpha.tsv", 
                      sep="\t", row.names=1, header=T)
#load in number of observed amplicons
num_ASVs <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/diversity/alpha_data/observed_asvs.tsv", 
                       sep="\t", row.names = 1, header=T)
#load in evenness scores
evenness <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/diversity/alpha_data/evenness.tsv", 
                       sep="\t", row.names = 1, header=T)
#check if all the row names from alpha div match up
rownames(faiths_pd) == rownames(shannon)
rownames(num_ASVs) == rownames(faiths_pd)
rownames(evenness) == rownames(faiths_pd)
#okay all are the same

#okay match up metadata data to faiths pd scores
metadata_fix <- metadata[rownames(faiths_pd),]

#attach alpha scores to the metadata DF
metadata_fix$faiths_pd <- faiths_pd$faith_pd
metadata_fix$shannon <- shannon$shannon
metadata_fix$num_ASVs <- num_ASVs$observed_otus
metadata_fix$evenness <- evenness$pielou_e

#check if they were input correctly
str(metadata_fix$faiths_pd)

#okay lets check distrubtion of data
hist(metadata_fix$faiths_pd, breaks=20)

#looks like a non-normal distrubition lets double check with a shapiro test
shapiro.test(metadata_fix$num_ASVs)
#distrubition is not normal will need to use non-parameteric tests

#test faiths_pd between IC and NIC
wilcox.test(faiths_pd ~ Infection_in_6, data=metadata_fix)
#looks to be significantly different

wilcox.test(shannon ~ Infection_in_6, data=metadata_fix)

#exact p-value can't be determined due to ties. (use exact = F to get a estimate of the p value)
wilcox.test(num_ASVs ~ Infection_in_6, data=metadata_fix, exact=F)
wilcox.test(evenness ~ Infection_in_6, data=metadata_fix)
#others are not significantly different, although num_ASVs is close... 

#okay lets check if ABX exposures are significantly different

#vanco exposure
wilcox.test(faiths_pd ~ Vanco_exposure, data=metadata_fix)
#anti-fungal exposure
wilcox.test(faiths_pd ~ anti_fungal_exposure, data=metadata_fix)

#piptaz exposure
wilcox.test(faiths_pd ~ PIP_TAZ_exposure, data=metadata_fix)
#other abx
wilcox.test(faiths_pd ~ OTHER_ABX, data=metadata_fix)

wilcox.test(faiths_pd ~ ABX_OT_PIPTaz, data=metadata_fix)

### AGE 
### problem is that if you do age in months this takes into account both the age at diagnosis
##but also the time throughout therapy...
cor.test(metadata_fix$faiths_pd, metadata_fix$Age_in_months)
plot(metadata_fix$faiths_pd ~ metadata_fix$Age_in_months)

cor.test(metadata_fix$faiths_pd, metadata_fix$Age_at_diagnosis)
plot(metadata_fix$faiths_pd ~ metadata_fix$Age_at_diagnosis)
library(ggplot2)
library(cowplot)
age_faiths <- ggplot(metadata_fix, aes(Age_at_diagnosis, faiths_pd, color=Infection_Col)) + geom_point(size=5) +
  xlab("Age at Diagnosis") + ylab("Faith's Phylogenetic Diversity") + guides(size=F, color=guide_legend(title="Infectious Complications")) +
  scale_color_identity(guide = T, labels=c("No", "Yes"))
age_faiths


### highly correlated with one another. 
cor.test(metadata_fix$Age_at_diagnosis, metadata_fix$faiths_pd, method="spearman")
cor.test(metadata_fix$Days_since_therapy_start, metadata_fix$faiths_pd, method="spearman")


#hmm this is a concerning figure as it could indicate our results...
age_months <- ggplot(metadata_fix, aes(Age_in_months, faiths_pd, color=Infection_Col)) + geom_point(size=5) +
  ylab("Faith's Phylogenetic Diversity") + xlab("Age in Months at Sample Collection") +guides(size=F, 
                                                                                              color=guide_legend(title="Infectious Complications")) +
  scale_color_identity(guide=T, labels=c("No", "Yes"))
age_months
#plot for faith's pd and abx exposure
library(ggplot2)

library(cowplot)
fill_colors_vanco <- ifelse(metadata_fix$Vanco_exposure=="N", "#FFF68F", "#800080")
faiths_plot_vanco <- ggplot(metadata_fix, aes(Vanco_exposure, faiths_pd, fill=fill_colors_vanco)) + geom_boxplot()+ 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Faith's Phylogenetic Diversity")+ 
  guides(fill=F, color=F)+ scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + scale_fill_identity() + xlab("Vancomycin Exposure")
faiths_plot_vanco

fill_colors_fung <- ifelse(metadata_fix$anti_fungal_exposure=="N", "#FFF68F", "#800080")
faiths_plot_fungal <- ggplot(metadata_fix, aes(anti_fungal_exposure, faiths_pd, fill=fill_colors_fung)) + geom_boxplot()+ 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Faith's Phylogenetic Diversity")+ 
  guides(fill=F, color=F)+ scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + scale_fill_identity() + xlab("Anti-fungal Exposure")
faiths_plot_fungal
boxplot(faiths_pd ~ anti_fungal_exposure, data=metadata_fix)

fill_colors_pip <- ifelse(metadata_fix$PIP_TAZ_exposure=="N", "#FFF68F", "#800080")
faiths_plot_piptaz <- ggplot(metadata_fix, aes(PIP_TAZ_exposure, faiths_pd, fill=fill_colors_pip)) + geom_boxplot()+ 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Faith's Phylogenetic Diversity")+ 
  guides(fill=F, color=F)+ scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + scale_fill_identity() + xlab("Piperacillin Tazobactam Exposure")
faiths_plot_piptaz

fill_colors_other <- ifelse(metadata_fix$OTHER_ABX=="N", "#FFF68F", "#800080")
faiths_plot_other_ABX <- ggplot(metadata_fix, aes(OTHER_ABX, faiths_pd, fill=fill_colors_other)) + geom_boxplot()+ 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Faith's Phylogenetic Diversity")+ 
  guides(fill=F, color=F)+ scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + scale_fill_identity() + xlab("Other Antibiotic Exposure")
faiths_plot_other_ABX

#combined plots

ABX_alpha_final_plot <- plot_grid(faiths_plot_vanco, faiths_plot_fungal, faiths_plot_piptaz, faiths_plot_other_ABX, labels = 'AUTO')
ABX_alpha_final_plot

#lets check how well infection and vanco exposure co-occur
fisher.test(metadata_fix$Infection_in_6, metadata_fix$Vanco_exposure)
#highly colrrelated with each other.....

#lets check how abx other then pip-taz is related to ABX
fisher.test(metadata_fix$Infection_in_6, metadata_fix$ABX_OT_PIPTaz)


#check for a significant difference overtime
cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$faiths_pd, method="spearman", exact=F)
#close to significant... (using peason's correlation)

#test other four diversity metrics
cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$shannon, method="spearman", exact=F)
cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$evenness, method="spearman", exact=F)
cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$num_ASVs, method="spearman", exact=F)
#okay lets make a logisitic model to try and classify samples from NIC and IC patients
#okay only significant correlation is shannon 

mod1 <- glm(Infection_in_6 ~ scale(faiths_pd), data=metadata_fix, family=binomial())
summary(mod1)
#faith's pd look to be a significant predictor

#try controlling for the time point of collection and treatment type of individual
mod2 <- glm(Infection_in_6 ~ scale(faiths_pd) + scale(Days_since_therapy_start) + Treatment_type, data=metadata_fix, family=binomial())
summary(mod2)
#model doesn't look to be significantly different from the first (faiths_pd is still the only significant predictor)

mod3 <- glm(Infection_in_6 ~ scale(faiths_pd) + scale(Days_since_therapy_start) + Treatment_type + Age_at_diagnosis, 
            data=metadata_fix, family=binomial())
mod3



test <- lm(faiths_pd ~ Infection_in_6 + Age_at_diagnosis + Days_since_therapy_start, data=metadata_fix)
summary(test)
summary(mod3)
anova(mod1, mod2, test="Chisq")
#not significantly different.... 

#generate a ROC for mod2
probs <- predict(mod3, type=c("response"))

library(pROC)
g <- roc(Infection_in_6 ~ probs, data=metadata_fix)
plot(g)
g
AUC <- 0.7676
#get kappa stat
library(caret)
threshold <- 0.5
confusionMatrix(factor(probs > threshold), factor(metadata_fix$Infection_in_6=="Y"), positive="TRUE")


#check for differences in faiths_pd due to sex (note all females in study are IC)
wilcox.test(faiths_pd ~ Sex, data=metadata_fix)
#not significantly different between the two groups

#### done data analysis on alpha differences between IC and NIC 

#make plots for this data:
library(ggplot2)
library(cowplot)

#make Figure1
faiths_plot <- ggplot(metadata_fix, aes(Infection_in_6, faiths_pd, fill=metadata_fix$Infection_Col)) + geom_boxplot()+ 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Faith's Phylogenetic Diversity")+ 
  guides(fill=F, color=F)+ scale_x_discrete(labels=c("No"="No", "Y"="Yes")) + scale_fill_identity()
faiths_plot

shannon_plot <- ggplot(metadata_fix, aes(Infection_in_6, shannon, fill=metadata_fix$Infection_Col)) + geom_boxplot() +
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Shannon Diversity") + 
  guides(fill=F, color=F) + scale_x_discrete(labels=c("No"="No", "Y"="Yes")) + scale_fill_identity()
shannon_plot

num_ASVs_plot <- ggplot(metadata_fix, aes(Infection_in_6, num_ASVs, fill=metadata_fix$Infection_Col)) + geom_boxplot() +
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Amplicon Sequence Variantes") + 
  guides(fill=F, color=F) + scale_x_discrete(labels=c("No"="No", "Y"="Yes")) + scale_fill_identity()
num_ASVs_plot

evenness_plot <- ggplot(metadata_fix, aes(Infection_in_6, evenness, fill=metadata_fix$Infection_Col)) + geom_boxplot() +
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Evenness") + scale_x_discrete(labels=c("No"="No", "Y"="Yes")) + 
  guides(fill=F, color=F) + scale_fill_identity()
evenness_plot

#generate figure1
#note that significant stars were added in post using libreoffice draw
figure1 <- plot_grid(faiths_plot, num_ASVs_plot, evenness_plot, shannon_plot, labels = "AUTO", nrow=2) + draw_label("Infectious Complications", y=0.03)
figure1


#generate Sup Fig 2
ROC <- plot(g)

library(gridGraphics)


#plot ROC curve

ggroc <- function(roc, showAUC = TRUE, interval = 0.2, breaks = seq(0, 1, interval)){
  require(pROC)
  if(class(roc) != "roc")
    simpleError("Please provide roc object from pROC package")
  plotx <- rev(roc$specificities)
  ploty <- rev(roc$sensitivities)
  
  ggplot(NULL, aes(x = plotx, y = ploty)) + 
    geom_step() +
    scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) + 
    scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    theme_bw() + 
    theme(axis.ticks = element_line(color = "grey80")) +
    coord_equal() + 
    annotate("text", x = interval/1.25, y = interval/2, vjust = 0, label = paste("AUC =",sprintf("%.2f",roc$auc)))
}
#put all the plots together

SupFig_2<- ggroc(g)
SupFig_2

### analysis pre vs post vs never data

#okay lets do pre v post
#subset data to look at pre and post
metadata_fix$Pre_Post_never_infect
detach()
metadata_pre_post <- metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Pre" | metadata_fix$Pre_Post_never_infect=="Post"),]
metadata_pre_post$Pre_Post_never_infect <- droplevels(metadata_pre_post$Pre_Post_never_infect)
attach(metadata_pre_post)

##### We could look at further into pre vs never vs post...
wilcox.test(faiths_pd ~ Pre_Post_never_infect, data=metadata_pre_post)
wilcox.test(shannon ~ Pre_Post_never_infect, data=metadata_pre_post)
wilcox.test(evenness ~ Pre_Post_never_infect, data=metadata_pre_post)
wilcox.test(num_ASVs ~ Pre_Post_never_infect, data=metadata_pre_post)
#okay we find nothing of significances.
detach()

#####################################################################
#lets look at pre vs never
metadata_pre_never <- metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Pre" | metadata_fix$Pre_Post_never_infect=="Never"),]
metadata_pre_never$Pre_Post_never_infect <- droplevels(metadata_pre_never$Pre_Post_never_infect)
wilcox.test(faiths_pd ~ Pre_Post_never_infect, data=metadata_pre_never)
wilcox.test(shannon ~ Pre_Post_never_infect, data=metadata_pre_never)
wilcox.test(evenness ~ Pre_Post_never_infect, data=metadata_pre_never)
wilcox.test(num_ASVs ~ Pre_Post_never_infect, data=metadata_pre_never, exact=F)
#okay we find nothing of significances

#test post vs never
metadata_post_never <- metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Post" | metadata_fix$Pre_Post_never_infect=="Never"),]
metadata_post_never$Pre_Post_never_infect <- droplevels(metadata_post_never$Pre_Post_never_infect)

wilcox.test(faiths_pd ~ Pre_Post_never_infect, data=metadata_post_never)
wilcox.test(shannon ~ Pre_Post_never_infect, data=metadata_post_never)
wilcox.test(num_ASVs ~ Pre_Post_never_infect, data=metadata_post_never, exact=F)
wilcox.test(evenness ~ Pre_Post_never_infect, data=metadata_post_never)

## faiths_pd significantly different and number of ASVs close to be significantly different.

post_never_faiths <- ggplot(metadata_fix, aes(Pre_Post_never_infect, faiths_pd, fill=metadata_fix$State_Cols)) + 
  geom_boxplot() + geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) +
  xlab("") + ylab("Faith's Phylogenetic Diversity")+ 
  guides(fill=F, color=F) + scale_fill_identity()
post_never_faiths


post_never_evenness <- ggplot(metadata_fix, aes(Pre_Post_never_infect, evenness, fill=metadata_fix$State_Cols)) + 
  geom_boxplot() + geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) +
  xlab("") + ylab("Evenness")+ 
  guides(fill=F, color=F) + scale_fill_identity()
post_never_evenness


post_never_shannon <- ggplot(metadata_fix, aes(Pre_Post_never_infect, shannon, fill=metadata_fix$State_Cols)) + 
  geom_boxplot() + geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) +
  xlab("") + ylab("Shannon Diversity")+ 
  guides(fill=F, color=F) + scale_fill_identity()
post_never_shannon


post_never_ASVs <-  ggplot(metadata_fix, aes(Pre_Post_never_infect, num_ASVs, fill=metadata_fix$State_Cols)) + 
  geom_boxplot() + geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) +
  xlab("") + ylab("Amplicon Sequence Variants")+ 
  guides(fill=F, color=F) + scale_fill_identity()
post_never_ASVs

Sup_fig_3 <- plot_grid(post_never_shannon, post_never_ASVs, post_never_evenness, post_never_faiths, labels = "AUTO")
Sup_fig_3


fisher.test(metadata_post_never$Pre_Post_never_infect, metadata_post_never$Sex)
#again significant were all females are post...
fisher.test(metadata_post_never$Pre_Post_never_infect, metadata_post_never$Vanco_exposure)
#against highly associated with each other
wilcox.test(faiths_pd ~ Vanco_exposure, data=metadata_post_never, exact=F)
#vanco exposure not significantly associated 
wilcox.test(faiths_pd ~ ABX_OT_PIPTaz, data=metadata_post_never, exact=F)
wilcox.test(faiths_pd ~ anti_fungal_exposure, data=metadata_post_never, exact=F)
wilcox.test(faiths_pd ~ PIP_TAZ_exposure, data=metadata_post_never, exact=F)
wilcox.test(faiths_pd ~ OTHER_ABX, data=metadata_post_never, exact=F)

#look at sex differences at baseline
metadata_base <- metadata_fix[which(metadata_fix$BASELINE_S=="Y"),]
wilcox.test(faiths_pd ~ Sex, data=metadata_base)

# analysis of differences over time.

cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$faiths_pd, method="pearson", exact=F)
#Significant

#test other four diversity metrics
cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$shannon, method="spearman", exact=F)
cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$evenness, method="spearman", exact=F)
cor.test(x=metadata_fix$Days_since_therapy_start, y=metadata_fix$num_ASVs, method="spearman", exact=F)
#close to significant.

#plot these.

#Could add regression line or cor line???
make_alpha_over_time <- function(x, alpha_metric, title){
    ret_plot <- ggplot(x, aes(Days_since_therapy_start, y=alpha_metric, fill=P_Color)) + geom_point(shape=21, size=3) +
      ylab(title) + xlab("Days Post Treatment Start") + guides(fill=F)
    return(ret_plot)
}

faith_OT_time <- make_alpha_over_time(metadata_fix, metadata_fix$faiths_pd, "Faith's Phylogenetic Diversity")
faith_OT_time
shannon_OT_time <- make_alpha_over_time(metadata_fix, metadata_fix$shannon, "Shannon Diversity")
evenness_OT_time <- make_alpha_over_time(metadata_fix, metadata_fix$evenness, "Evenness")
num_ASVs_OT_time <- make_alpha_over_time(metadata_fix, metadata_fix$num_ASVs, "Observed ASVs")

alpha_OT_time_Comb <- plot_grid(faith_OT_time, shannon_OT_time, evenness_OT_time, evenness_OT_time, num_ASVs_OT_time, labels="AUTO")
alpha_OT_time_Comb


#### Individuals patients overtime.


#P10 overtime

metadata_P10 <- metadata_fix[which(metadata_fix$Patient=='P10'),]

P10_OT_time <- ggplot(metadata_P10, aes(Days_since_therapy_start, y=faiths_pd, fill=P_Color)) + geom_point(size=3, shape=21) +
  ylab("") +xlab("Days Post Treatment Start") + scale_fill_identity() + guides(fill=F) + ggtitle("P10")
P10_OT_time

metadata_P13 <- metadata_fix[which(metadata_fix$Patient=="P13"),]

P13_OT_time <- ggplot(metadata_P13, aes(Days_since_therapy_start, y=faiths_pd, fill=P_Color)) + geom_point(size=3, shape=21) +
  ylab("Faith's Phylogenetic Diversity") + xlab("Days Post Treatment Start") + scale_fill_identity() + guides(fill=F) + ggtitle("P13")
P13_OT_time

metadata_P16 <- metadata_fix[which(metadata_fix$Patient=="P16"),]

P16_OT_time <- ggplot(metadata_P16, aes(Days_since_therapy_start, y=faiths_pd, fill=P_Color)) + geom_point(size=3, shape=21) +
  ylab("") + xlab("Days Post Treatment Start") + scale_fill_identity() + guides(fill=F) + ggtitle("P16")
P16_OT_time


Case_PD_OT_time <- plot_grid(P10_OT_time, P13_OT_time, P16_OT_time, labels = "AUTO", cols = 1)
Case_PD_OT_time


metadata_P17 <- metadata_fix[which(metadata_fix$Patient=="P17"),]


#most of P17 samples didn't sequence well so can't do
P17_OT_time <- ggplot(metadata_P17, aes(Days_since_therapy_start, y=faiths_pd, fill=P_Color)) + geom_point(size=3, shape=21) +
  ylab("Faith's Phylogenetic Diversity") + xlab("Days Post Treatment Start") + scale_fill_identity() + guides(fill=F)
P17_OT_time



#### Baseline analysis for MRD

metadata_BASELINE <- metadata_fix[which(metadata_fix$BASELINE_S=="Y"),]
wilcox.test(metadata_BASELINE$faiths_pd ~ metadata_BASELINE$Infection_in_6)
boxplot(metadata_BASELINE$faiths_pd ~ metadata_BASELINE$Infection_in_6)

wilcox.test(metadata_BASELINE$faiths_pd ~ metadata_BASELINE$MRD)

table(metadata_BASELINE$Infection_in_6, metadata_BASELINE$PIP_TAZ_exposure)
metadata_BASELINE$Infection_in_6

### test subsequent infection's


metadata_sub <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17_subsequent.csv", sep="\t",
                       row.names=1, header=T, comment.char="")
metadata_sub_fix <- droplevels(metadata_sub[rownames(metadata_fix),])
rownames(metadata_sub_fix) == rownames(metadata_fix)

metadata_sub_fix$faiths_pd <- metadata_fix$faiths_pd

wilcox.test(metadata_sub_fix$faiths_pd ~ metadata_sub_fix$Subsequent_GI)
wilcox.test(metadata_sub_fix$faiths_pd ~ metadata_sub_fix$Subsequent_BSI)

#make fig
fill_colors_sub_BSI <- ifelse(metadata_sub_fix$Subsequent_BSI=="No", "#FFF68F", "#800080")
sub_BSI_div <- ggplot(metadata_sub_fix, aes(Subsequent_BSI, faiths_pd, fill=fill_colors_sub_BSI)) + geom_boxplot()+ 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_fix$P_Color)) + xlab("") + ylab("Faith's Phylogenetic Diversity")+ 
  guides(fill=F, color=F)+ scale_x_discrete(labels=c("N"="No", "Y"="Yes")) + scale_fill_identity() + xlab("Subsequent Blood Stream Infection")
sub_BSI_div

wilcox.test(metadata_sub_fix$faiths_pd ~ metadata_sub_fix$Subsequent_FN)
wilcox.test(metadata_sub_fix$faiths_pd ~ metadata_sub_fix$Infection_in_6)

metadata_sub_baseline <- metadata_sub[rownames(metadata_BASELINE),]
metadata_sub_baseline$faiths_pd <- metadata_BASELINE$faiths_pd

wilcox.test(metadata_sub_baseline$faiths_pd ~ metadata_sub_baseline$Subsequent_BSI)
wilcox.test(metadata_sub_baseline$faiths_pd ~ metadata_sub_baseline$Subsequent_GI)
wilcox.test(metadata_sub_baseline$faiths_pd ~ metadata_sub_baseline$Subsequent_FN)






metadata_sub_baseline$Subsequent_BSI
boxplot(metadata_sub_fix$faiths_pd ~ metadata_sub_fix$Subsequent_BSI)

