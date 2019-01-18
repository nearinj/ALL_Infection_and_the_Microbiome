#### Do Hakim Analysis
### Baseline analysis based on increased level of Proteobacteria
### Does increased proteo at baseline mean increased risk for infection
## load in table

taxa <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/phylum_table.tsv",
                   sep="\t", skip=1, comment.char="", quote="", header=T, row.names = 1)

#read in metadata to get baseline samples and find out infection vs no infection

#convert to RA

taxa <- sweep(taxa, 2, colSums(taxa), '/')

#read in baseline

metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17_subsequent.csv",
                       sep="\t", comment.char = "", quote = "", header=T, row.names = 1)

rownames(metadata) <- gsub("-", ".", rownames(metadata))


#only keep baseline samples
metadata_filt <- metadata[colnames(taxa),]
#get baseline samples
metadata_baseline <- metadata_filt[which(metadata_filt$BASELINE_S=="Y"),]

#11 baseline samples


#okay now clean up taxa table

taxa_baseline <- taxa[,rownames(metadata_baseline)]
taxa_baseline <- data.frame(t(taxa_baseline))

boxplot(taxa_baseline$k__Bacteria.p__Proteobacteria ~ metadata_baseline$Infection_in_6)
metadata_baseline$Infection_in_6 <- ifelse(metadata_baseline$Infection_in_6=="Y", "Yes", "No")
#significantly different at baseline interesting....
wilcox.test(taxa_baseline$k__Bacteria.p__Proteobacteria ~ metadata_baseline$Infection_in_6)
library(cowplot)
library(ggplot2)
baseline_proteo <- ggplot(metadata_baseline, aes(metadata_baseline$Infection_in_6, y=taxa_baseline$k__Bacteria.p__Proteobacteria, 
                                                 fill=Infection_Col)) + geom_boxplot() + scale_fill_identity() +
  geom_point(aes(fill=P_Color), size=5, shape=21) +xlab("Infectious Complications") + ylab("Relative Abudance of Proteobacteria") +
  scale_y_continuous(labels = scales::percent)

baseline_proteo
#okay we need to do a grey test now... So we need the number of days until each sample had their first infection for these baseline samples.


##Baseline graph.


#okay lets do a grey test

library(cmprsk)

metadata_baseline$Infection_in_6

summary(metadata_baseline$Infection_in_6)
#okay in their problem they grouped by a RA greater than 0.01%

#divide grp by proteobacteria

#0.01% doesn't work .... everyone is placed as high proteobacteria levels we will adjust this
#we set the level to 1% and it works much better
taxa_baseline$pro_grp <- ifelse(taxa_baseline$k__Bacteria.p__Proteobacteria >= 0.01, "High", "Low")

taxa_baseline$pro_grp <- factor(taxa_baseline$pro_grp)
summary(taxa_baseline$pro_grp, metadata_baseline$Infection_in_6)

fit=cuminc(metadata_baseline$Time_to_infection_Base, metadata_baseline$Infection_in_6, taxa_baseline$pro_grp, cencode = "N")
rownames(metadata_baseline) == rownames(taxa_baseline)

fit
plot(fit, xlab="Days")

test <- cbind(taxa_baseline, metadata_baseline)
table(test$pro_grp, test$Infection_in_6)
#okay so all invdidivuals that had proteobacteria under 0.01% didn't have infections
#5/7 individuals with RA above 0.1% had infection
#do to low number of patients we cannot look at gray test....


#lets see if we get similar results from metagenomic data

MGS_phylum <- read.table("~/projects/MALL-MGS-ALL/Phylum_metaphlan2.tsv", sep="\t", header=T, quote="", row.names = 1)
metadata_baseline_MGs <- metadata[which(metadata$BASELINE_S=="Y"),]
MGS_phylum_base <- MGS_phylum[,rownames(metadata_baseline_MGs)]
MGS_phylum_base <- data.frame(t(MGS_phylum_base))
#okay first lets do a wilcoxon test

wilcox.test(MGS_phylum_base$k__Bacteria.p__Proteobacteria ~ metadata_baseline_MGs$Infection_in_6)
#RA of proteobacteria is significantly different.

#okay divide into high and low and see what it looks like.
MGS_phylum_base$pro_grp <- ifelse(MGS_phylum_base$k__Bacteria.p__Proteobacteria>=1, "High", "Low")
#0.01% cut off doesn't work in MGS
MGS_phylum_base$pro_grp <- factor(MGS_phylum_base$pro_grp)
rownames(MGS_phylum_base) == rownames(metadata_baseline_MGs)
test_mgs <- cbind(MGS_phylum_base, metadata_baseline_MGs)

table(test_mgs$pro_grp, test_mgs$Infection_in_6)

#very similar results.

#okay lets make tables for these results... and add to results section of our manuscript!


fisher.test(test_mgs$Infection_in_6, test_mgs$pro_grp)

test_table <- glm(Infection_in_6 ~ k__Bacteria.p__Proteobacteria, data=test_mgs, family=binomial())
test_table
summary(test_table)

library(pROC)

probs <- predict(test_table, type=c("response"))
probs
