#### Analysis of Metaphlan Output

### load in metadata
metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17.csv", sep="\t", 
                       header=T, row.names = 1, comment.char="")
rownames(metadata) <- gsub("-",".",rownames(metadata))

#only get metadata from samples within 6 Months
metadata_fix <- metadata[which(metadata$Within_183_days=="Y"),]

### load in speices table
spec_table <- read.table("~/projects/MALL-MGS-ALL/species_metaphlan2.tsv", sep="\t", header=T, row.names = 1)
#flip table
spec_table_flip <- as.data.frame(t(spec_table))

#match metadata rownames to species table 
spec_table_fix <- spec_table_flip[rownames(metadata_fix),]

#check for significantly different species
#run wilcoxon across all taxa 
pvals <- apply(spec_table_fix, 2, function(x) wilcox.test(x ~ metadata_fix$Infection_in_6, exact=F)$p.value)
#correct for multiple tests
qval <- p.adjust(pvals, "fdr")
#get those that meet our alpha value criteria
hits <- qval < 0.05
which(hits)

#load in table contain genera
gen_tab <- read.table("~/projects/MALL-MGS-ALL/test/Genera_table.txt", sep="\t", header=T, row.names = 1)
gen_tab_flip <- as.data.frame(t(gen_tab))
gen_tab_fix <- gen_tab_flip[rownames(metadata_fix),]

pvals_gen <- apply(gen_tab_fix, 2, function(x) wilcox.test(x ~ metadata_fix$Infection_in_6, exact=F)$p.value)
qval_gen <- p.adjust(pvals_gen, "fdr")
hits_gen <- qval_gen < 0.05
which(hits_gen)





#6 significant hits from species table.
#pull them out into a new table
sig_specs <- spec_table_fix[,which(hits)]
#calculate means of infection and non-infection

Infection_means <- colMeans(sig_specs[which(metadata_fix$Infection_in_6=="Y"),])
Infection_means

Non_Infection_means <- colMeans(sig_specs[which(metadata_fix$Infection_in_6=="No"),])
Non_Infection_means

Infection_means - Non_Infection_means
#bind metadata to the signicant species table
metadata_sig_spec <- cbind(metadata_fix, sig_specs)

metadata_sig_spec$Infection_in_6 <- ifelse(metadata_sig_spec$Infection_in_6=="Y", "Yes", "No")



library(ggplot2)
#make Figure 3
panA <- ggplot(metadata_sig_spec, aes(Infection_in_6, sig_specs[,1], fill=metadata_sig_spec$Infection_Col)) + 
  geom_boxplot(outlier.shape = NA) +
  xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_sig_spec$P_Color)) + 
  ggtitle(expression(paste(italic('Fecalibacterium prausnitzii')))) +
  scale_fill_identity()
panA

panB <- ggplot(metadata_sig_spec, aes(Infection_in_6, sig_specs[,2], fill=metadata_sig_spec$Infection_Col)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_sig_spec$P_Color)) +
  ggtitle(expression(paste(italic('Brevundimonas diminuta')))) +
  scale_fill_identity()

panB

panC <- ggplot(metadata_sig_spec, aes(Infection_in_6, sig_specs[,3], fill=metadata_sig_spec$Infection_Col)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_sig_spec$P_Color)) + 
  ggtitle(expression(paste(italic('Agrobacterium tumefaciens')))) +
  scale_fill_identity()
panC

panD <- ggplot(metadata_sig_spec, aes(Infection_in_6, sig_specs[,4], fill=metadata_sig_spec$Infection_Col)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_sig_spec$P_Color)) +
  ggtitle(expression(paste(italic("Agrobacterium "),"unclassified"))) +
  scale_fill_identity()
panD

#make a 3 by 2 boxplot panel for this figure
panE <- ggplot(metadata_sig_spec, aes(Infection_in_6, sig_specs[,5], fill=metadata_sig_spec$Infection_Col)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_sig_spec$P_Color)) +
  ggtitle(expression(paste(italic("Achromobacter")," unclassified"))) +
  scale_fill_identity()
panE

panF <- ggplot(metadata_sig_spec, aes(Infection_in_6, sig_specs[,6], fill=metadata_sig_spec$Infection_Col)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_sig_spec$P_Color)) +
  ggtitle(expression(paste(italic("Alcaligenes "), "unclassified"))) +
  scale_fill_identity()
panF


library(cowplot)

Figure3 <- plot_grid(panA, panB, panC, panD, panE, panF, nrow=2, ncol=3, labels="AUTO")
Figure3


#generate genera figure (sup figure 5)

sig_genera <- gen_tab_fix[,which(hits_gen)]
colnames(sig_genera) <- gsub(".*g__","",colnames(sig_genera))
metadata_sig_gen <- cbind(metadata_fix, sig_genera)

make_boxplots <- function(x){
  plot <- ggplot(metadata_sig_gen, aes(Infection_in_6, sig_genera[,x], fill=metadata_sig_gen$Infection_Col))+ geom_boxplot(outlier.shape = NA) +
    xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
    geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_sig_gen$P_Color)) +
    ggtitle(colnames(sig_genera)[x]) +
    scale_fill_identity()
  return(plot)
}
plot_list_genera <- list()
for (i in 1:ncol(sig_genera)){
  plot_list_genera[[i]] <- make_boxplots(i)
}

Sup_fig_5 <- plot_grid(plotlist = plot_list_genera, labels="AUTO")  
Sup_fig_5
plot_list_genera[[2]]
#okay done these plots.


#Do pre post never analysis

metadata_pre_post <- droplevels(metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Pre" | metadata_fix$Pre_Post_never_infect=="Post"),])
spec_tab_pre_post <- spec_table_fix[rownames(metadata_pre_post),]

pval_pre_post <- apply(spec_tab_pre_post, 2, function(x) wilcox.test(x ~ metadata_pre_post$Pre_Post_never_infect, exact=F)$p.value)
qval_pre_post <- p.adjust(pval_pre_post, "fdr")
hit_pre_post <- qval_pre_post < 0.05
which(hit_pre_post)
#no significant hits

metadata_pre_never <- droplevels(metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Pre" | metadata_fix$Pre_Post_never_infect=="Never"),])
spec_tab_pre_never <- spec_table_fix[rownames(metadata_pre_never),]
pval_pre_never <- apply(spec_tab_pre_never, 2, function(x) wilcox.test(x ~ metadata_pre_never$Pre_Post_never_infect, exact=F)$p.value)
qval_pre_never <- p.adjust(pval_pre_never, "fdr")
hit_pre_never <- qval_pre_never < 0.05
which(hit_pre_never)
#no significant hits

metadata_post_never <- droplevels(metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Post" | metadata_fix$Pre_Post_never_infect=="Never"),])
spec_tab_post_never <- spec_table_fix[rownames(metadata_post_never),]
pval_tab_post_never <- apply(spec_tab_post_never, 2, function(x) wilcox.test(x ~ metadata_post_never$Pre_Post_never_infect, exact=F)$p.value)
qval_tab_post_never <- p.adjust(pval_tab_post_never, "fdr")
hit_post_never <- qval_tab_post_never < 0.05
which(hit_post_never)

# we get the same hits but with an extra... this is most likely due to less fdr correction...
# so what we really are seeing in difference is never vs infection we don't really see differences in pre
# pre seems to hover around in between the two.



#post_never plots
sig_spec_post_never <- spec_tab_post_never[,which(hit_post_never)]

panA_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,1], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Fecalibacterium prausnitzii"))))
panA_post_never

panB_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,2], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Brevundimonas diminuta"))))
panB_post_never

panC_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,3], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Agrobacterium tumefaciens"))))
panC_post_never

panD_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,4], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Agrobacterium "),"unclassified")))
panD_post_never

panG_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,5], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Achromobacter piechaudii"))))
panG_post_never

panE_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,6], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("Relative Abundance") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Achromobacter "), "unclassified")))
panE_post_never

panF_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,7], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("RA of Alcaligenes unclassified") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Alcaligenes "),"unclassified")))
panF_post_never

panH_post_never <- ggplot(metadata_post_never, aes(Pre_Post_never_infect, sig_spec_post_never[,8], fill=metadata_post_never$State_Cols)) + geom_boxplot(outlier.shape = NA) +
  xlab("Infectious State") + ylab("RA of Bordetella unclassified") + guides(fill=F) + 
  geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=metadata_post_never$P_Color)) +theme(axis.title = element_text(size=12)) +
  scale_fill_identity() +
  ggtitle(expression(paste(italic("Bordetella "),"unclassified")))
panH_post_never

#post never panel (Sup Fig 6)
post_never_panel <- plot_grid(panA_post_never,
                              panB_post_never,
                              panC_post_never,
                              panE_post_never,
                              panF_post_never,
                              panG_post_never,
                              panH_post_never)
post_never_panel

