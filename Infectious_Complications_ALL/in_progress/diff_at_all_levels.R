### script to find Differential abundant taxa from metaphlan tables

#read in species table
metaphlan_spec <- read.table('~/projects/MALL-MGS-ALL/species_metaphlan2.tsv',
                             sep="\t", header=T, comment.char="", quote="", row.names = 1)
#load in metadata table

metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17_subsequent.csv",
                       sep="\t", comment.char="", quote="", row.names = 1, header=T)
rownames(metadata) <- gsub("-",".",rownames(metadata))

metadata_fix <- metadata[which(metadata$Within_183_days=="Y"),]
metadata_fix$Infection_in_6 <- ifelse(metadata_fix$Infection_in_6=="Y","Yes","No")

#function to get differentailly abundant taxa

get_diff_taxa <- function(taxa_tab, metadata_info){
  
  taxa_tab <- data.frame(t(taxa_tab))
  
  taxa_tab_filt <- taxa_tab[rownames(metadata_info),]
  print(length(rownames(taxa_tab_filt)))
  print(length(colnames(taxa_tab_filt)))
  if(length(which(rownames(taxa_tab_filt) == rownames(metadata_info))) != length(rownames(taxa_tab_filt))){
    print("failure")
    return()
  }
  #okay run wilcoxon on each species
  
  pval <- apply(taxa_tab_filt, 2, function(x) wilcox.test(x ~ metadata_info$Infection_in_6, exact=F)$p.value)
  qval <- p.adjust(pval, "fdr")
  #make table of only differential abundant taxa
  diff_taxa <- taxa_tab_filt[,which(qval < 0.05),drop=F]

  return(diff_taxa)
}

Diff_species <- get_diff_taxa(metaphlan_spec, metadata_fix)



#load in other taxonomic levels for MGS data
#phylum level
Phylum_tab <- read.table("~/projects/MALL-MGS-ALL/Phylum_metaphlan2.tsv", 
                         sep="\t", comment.char="", quote="", row.names = 1, header=T)
colSums(Phylum_tab)
#class table
class_tab <- read.table("~/projects/MALL-MGS-ALL/class_metaphlan2.tsv",
                        sep="\t", comment.char="", quote="", row.names = 1, header=T)
colSums(class_tab)
#order table
order_tab <- read.table("~/projects/MALL-MGS-ALL/order_metaphlan2.tsv",
                        sep="\t", comment.char="", quote="", row.names = 1, header=T)
colSums(order_tab)
#family level
family_tab <- read.table("~/projects/MALL-MGS-ALL/family_metaphlan2.tsv",
                         sep="\t", comment.char="", quote="", row.names = 1, header=T)
colSums(family_tab)
#genus level

genus_tab <- read.table("~/projects/MALL-MGS-ALL/genus_metaphlan2.tsv", 
                        sep="\t", comment.char="", quote="", row.names =1, header=T)
colSums(genus_tab)
#okay now run diff on each

diff_phylum <- get_diff_taxa(Phylum_tab, metadata_fix)
colnames(diff_phylum)
diff_class <- get_diff_taxa(class_tab, metadata_fix)
diff_order <- get_diff_taxa(order_tab, metadata_fix)
diff_family <- get_diff_taxa(family_tab, metadata_fix)
diff_genus <- get_diff_taxa(genus_tab, metadata_fix)


#okay now need to write function that will make boxplots for each differential abundant taxa at each level
get_boxplots <- function(diff_taxa, metadata_fix, level){
  if(level=="g"){
    #drop any unclassified bacteria
    #diff_taxa <- diff_taxa[, -grep("__$", colnames(diff_taxa))]
    colnames(diff_taxa) <- gsub(".*g__","",colnames(diff_taxa))
  }else if(level=="s"){
    #drop any unclassified bacteria
    #diff_taxa <- diff_taxa[, -grep("__$", colnames(diff_taxa))]
    colnames(diff_taxa) <- gsub(".*s__", "", colnames(diff_taxa))
  } else if(level=="f"){
    #drop any unclassified bacteria
      #diff_taxa <- diff_taxa[, -grep("__$", colnames(diff_taxa))]
      colnames(diff_taxa) <- gsub(".*f__", "", colnames(diff_taxa))
  } else if(level=="o"){
    #drop any unclassified bacteria
      #diff_taxa <- diff_taxa[, -grep("__$", colnames(diff_taxa))]
      colnames(diff_taxa) <- gsub(".*o__", "", colnames(diff_taxa))
  } else if(level=="c"){
    #drop any unclassified bacteria
      #diff_taxa <- diff_taxa[, -grep("__$", colnames(diff_taxa))]
      colnames(diff_taxa) <- gsub(".*c__", "", colnames(diff_taxa))
  } else if(level=="p"){
      colnames(diff_taxa) <- gsub(".*p__","", colnames(diff_taxa))
  }
  
  colnames(diff_taxa) <- gsub("_"," ", colnames(diff_taxa))
  
  print(colnames(diff_taxa))
  plot_data <- cbind(metadata_fix, diff_taxa)
 
  plots <- list()
  
  for(i in 1:ncol(diff_taxa)){
    plots[[i]] <- make_boxplots(plot_data, diff_taxa, i)
  }
  
  return(plots)
}
#function that takes plot_data and makes ggplot out of it.
make_boxplots <- function(plot_data, diff_taxa, index){
  plot <- ggplot(plot_data, aes(Infection_in_6, diff_taxa[,index], fill=plot_data$Infection_Col))+ geom_boxplot(outlier.shape = NA) +
    xlab("Infectious Complications") + ylab("Relative Abundance") + guides(fill=F) + 
    geom_dotplot(stackdir = "center", binaxis = "y", aes(fill=plot_data$P_Color)) +
    ggtitle(colnames(diff_taxa)[index]) +
    scale_fill_identity()
  return(plot)
}

library(ggplot2)

library(cowplot)
Phylum_plots <- get_boxplots(diff_phylum, metadata_fix, "p")
print(Phylum_plots[[1]])
Phylum_plots <- plot_grid(plotlist=Phylum_plots, labels="AUTO")
Phylum_plots


#class plots
Class_plots <- get_boxplots(diff_class, metadata_fix, "c")
Class_plots <- plot_grid(plotlist=Class_plots, labels="AUTO")
Class_plots

#order plots
order_plots <- get_boxplots(diff_order, metadata_fix, "o")
order_plots <- plot_grid(plotlist=order_plots, labels="AUTO")
order_plots

#family plots
family_plots <- get_boxplots(diff_family, metadata_fix, "f")
family_plots <- plot_grid(plotlist=family_plots, labels="AUTO")
family_plots

#genus plots
genus_plots <- get_boxplots(diff_genus, metadata_fix, "g")
genus_plots <- plot_grid(plotlist=genus_plots, labels="AUTO")
genus_plots



#do pre vs never
metadata_pre_never <- metadata_fix[which(metadata_fix$Pre_Post_never_infect=="Pre" | metadata_fix$Pre_Post_never_infect=="Never"),]
#all done with all of these plots! not sure what level to exactly focus on here though....
pre_never_phylum <- get_diff_taxa(Phylum_tab, metadata_pre_never)
pre_never_phylum

pre_never_class <- get_diff_taxa(class_tab, metadata_pre_never)
colnames(pre_never_class)

pre_never_class_plot <- get_boxplots(pre_never_class, metadata_pre_never, "c")
pre_never_class_plot

pre_never_order <- get_diff_taxa(order_tab, metadata_pre_never)
colnames(pre_never_order)

pre_never_order_plot <- get_boxplots(pre_never_order, metadata_pre_never, "o")
pre_never_order_plot




pre_never_family <- get_diff_taxa(family_tab, metadata_pre_never)
colnames(pre_never_family)

pre_never_genus <- get_diff_taxa(genus_tab, metadata_pre_never)
colnames(pre_never_genus)
#do same analysis for 16S analysis
### load in 16S taxa tables

spec_tab_16 <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/spec_table.tsv",
                       sep="\t", skip=1, comment.char="", quote="", header=T, row.names=1)
colnames(spec_tab_16) <- gsub("-",".",colnames(spec_tab_16))

#convert to RA

spec_tab_16_RA <- sweep(spec_tab_16, 2, colSums(spec_tab_16), '/')*100
colSums(spec_tab_16_RA)

genus_tab_16 <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/genus_table.tsv",
                       sep="\t", skip=1, comment.char="", quote="", header=T, row.names=1)
colnames(genus_tab_16) <- gsub("-", ".", colnames(genus_tab_16))

#convert to RA

genus_tab_16_RA <- sweep(genus_tab_16, 2, colSums(genus_tab_16), '/')*100

family_tab_16 <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/family_table.tsv",
                       sep="\t", skip=1, comment.char="", quote="", header=T, row.names=1)
colnames(family_tab_16) <- gsub("-", ".", colnames(family_tab_16))

family_tab_16_RA <- sweep(family_tab_16, 2, colSums(family_tab_16), '/')*100


order_tab_16 <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/order_table.tsv",
                       sep="\t", skip=1, comment.char="", quote="", header=T, row.names=1)
colnames(order_tab_16) <- gsub("-", ".", colnames(order_tab_16))

order_tab_16_RA <- sweep(order_tab_16, 2, colSums(order_tab_16), '/')*100

class_tab_16 <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/class_table.tsv",
                       sep="\t", skip=1, comment.char="", quote="", header=T, row.names=1)
colnames(class_tab_16) <- gsub("-", ".", colnames(class_tab_16))

class_tab_16_RA <- sweep(class_tab_16, 2, colSums(class_tab_16), '/')*100

phylum_tab_16 <- read.table("~/projects/MALL-deblur/KulkarniMALL16S/6_mnth_updated_July/taxa/rare_taxa_tables/phylum_table.tsv",
                       sep="\t", skip=1, comment.char="", quote="", header=T, row.names=1)
colnames(phylum_tab_16) <- gsub("-", ".", colnames(phylum_tab_16))

phylum_tab_16_RA <- sweep(phylum_tab_16, 2, colSums(phylum_tab_16), '/')*100

metadata_fix_16 <- metadata[colnames(spec_tab_16),]

diff_spec_16 <- get_diff_taxa(spec_tab_16_RA, metadata_fix_16)

diff_genus_16 <- get_diff_taxa(genus_tab_16_RA, metadata_fix_16)

diff_class_16 <- get_diff_taxa(class_tab_16_RA, metadata_fix_16)

diff_phylum_16 <- get_diff_taxa(phylum_tab_16_RA, metadata_fix_16)

diff_order_16 <- get_diff_taxa(order_tab_16_RA, metadata_fix_16)

diff_family_16 <- get_diff_taxa(family_tab_16_RA, metadata_fix_16)


#make plots for each
spec_plot_16 <- get_boxplots(diff_spec_16, metadata_fix_16, "s")
spec_plot_16 <- plot_grid(plotlist=spec_plot_16, labels="AUTO")
spec_plot_16
#need to fix labels at top of plots.

genus_plot_16 <- get_boxplots(diff_genus_16, metadata_fix_16, "g")
genus_plot_16 <- plot_grid(plotlist=genus_plot_16, labels="AUTO")
genus_plot_16

family_plot_16 <- get_boxplots(diff_family_16, metadata_fix_16, "f")
family_plot_16 <- plot_grid(plotlist=family_plot_16, labels="AUTO")
family_plot_16

order_plot_16 <- get_boxplots(diff_order_16, metadata_fix_16, "o")
order_plot_16 <- plot_grid(plotlist=order_plot_16, labels="AUTO")
order_plot_16

class_plot_16 <- get_boxplots(diff_class_16, metadata_fix_16, "c")
class_plot_16 <- plot_grid(plotlist=class_plot_16, labels="AUTO")
class_plot_16

phylum_plot_16 <- get_boxplots(diff_phylum_16, metadata_fix_16, "p")
phylum_plot_16 <- plot_grid(plotlist=phylum_plot_16, labels="AUTO")
phylum_plot_16

### done

##### need to fix labels on all the plots and then determine how to summarize the data.... perhaps compare and contract the 16S and MGS data a bit???
