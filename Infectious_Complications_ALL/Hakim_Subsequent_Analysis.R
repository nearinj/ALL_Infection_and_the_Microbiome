#### Hakim Subsequent analysis


## load in family table


family_table <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/family_table.tsv",
                           skip=1, sep="\t", comment.char="", quote="", header=T, row.names = 1)

# convert to RA

family_table <- sweep(family_table, 2, colSums(family_table), '/')
colSums(family_table)
#converted to RA

#okay now read in metadata.

metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17_subsequent.csv",
                       sep="\t", header=T, quote="", comment.char = "", row.names = 1)
rownames(metadata) <- gsub("-", ".", rownames(metadata))


metadata_filt <- metadata[colnames(family_table),]

#okay calulate those with dominance by Enterococcaceae

family_table <- data.frame(t(family_table))
rowSums(family_table)

#no sample has 30% Enterococcaceae.... 
family_table$Entero_DOM <- ifelse(family_table$k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Enterococcaceae >= 0.30, 
                                  "Yes","No")
family_table$Entero_DOM <- factor(family_table$Entero_DOM)
summary(family_table$Entero_DOM)
#very few samples had family level assignment to enterococcaceae #i wonder if this has to do with primer bias

#lets look at Streptococcaceae 

family_table$Strepto_DOM <- ifelse(family_table$k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae >= .30, 
                                   "Yes", "No")
family_table$Strepto_DOM <- factor(family_table$Strepto_DOM)
summary(family_table$Strepto_DOM)
rownames(family_table) == rownames(metadata_filt)
test <- cbind(family_table, metadata_filt)
table(test$Strepto_DOM, test$Subsequent_GI)

fisher.test(test$Strepto_DOM, test$Subsequent_GI)
# 7 samples that are dom by Strepto
#calculate hazard ratio for subsequent GI, FN, and BSI
#do i need to go back and write the number of days until each event for each sample? #I don't think it makes sense here to include samples
#from the same individul does it?

#simply odds ratio test to see if that makes sense not significant... but this could be due to low sample numbers
fisher.test(family_table$Strepto_DOM, metadata_filt$Infection_in_6)

metadata_filt$Subsequent_GI<- droplevels(metadata_filt$Subsequent_GI)
fisher.test(family_table$Strepto_DOM, metadata_filt$Subsequent_GI)

metadata_filt$Subsequent_BSI <- droplevels(metadata_filt$Subsequent_BSI)
fisher.test(family_table$Strepto_DOM, metadata_filt$Subsequent_BSI)

metadata_filt$Subsequent_FN <- droplevels(metadata_filt$Subsequent_FN)
fisher.test(family_table$Strepto_DOM, metadata_filt$Subsequent_FN)

#none of them were significant. (confidence intervals are extreme between them)

#issues here.... #some samples may have been taken while they were already experiencing some sort of illness... 
#don't have the sample size to
#filter them out though...

wilcox.test(family_table$k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae ~ metadata_filt$Infection_in_6)
boxplot(family_table$k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae ~ metadata_filt$Infection_in_6)
#not sure how much would change if we took into account time using a hazard model such as the Anderson-Gill model...