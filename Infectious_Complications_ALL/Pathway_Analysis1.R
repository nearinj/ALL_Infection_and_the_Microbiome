### Pathway Analysis 


#load in metadata
metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17.csv", sep="\t", 
                       header=T, row.names=1, comment.char="")

rownames(metadata) <- gsub("-",".",rownames(metadata))

metadata_fix <- metadata[which(metadata$Within_183_days=="Y"),]

#load in pathway data
pathways <- read.table("~/projects/MALL-MGS-ALL/humann2_final_out/humann2_pathabundance_relab_unstratified.tsv", sep="\t", header=T, row.names=1, 
                       comment.char = "", quote="")
pathways_flip <- as.data.frame(t(pathways))
pathways_fix <- pathways_flip[rownames(metadata_fix),]

#run wilcoxon on each pathway
pvals <- apply(pathways_fix, 2, function(x) wilcox.test(x ~ metadata_fix$Infection_in_6, exact=F)$p.value)
#correct for flase discovery
qvals <- p.adjust(pvals, "fdr")
hits <- qvals < 0.05
which(hits)

#make DF contianing only the hits
sig_pathways <- pathways_fix[,which(hits)]
#52 significant pathways ... need to filter out human, and synthetic pathways

rownames(sig_pathways) <- gsub("_Abundance","",rownames(sig_pathways))
colnames(sig_pathways)

#add pathway types manually from metacyc
# 1 <- Interceronversions
# 2 <- Biolum
# 3 <- Biosynth
# 4 <- Degradtion
# 5 <- Detox
# 6 <- Precurser metabolites and energy
# 7 <- glycan pathways
# 8 <- macro mod
# 9 <- metabolic cluster
# 10 <- signal transduction
# 11 <- super pathways
# 99 <- human or synthetic pathway (removed in downstream steps)

sig_pathway_groupings <- c(1, 
                           4,
                           3,
                           3,
                           3,
                           3,
                           3,
                           4,
                           4,
                           4,
                           4,
                           6,
                           3,
                           3,
                           3,
                           3,
                           3,
                           6,
                           4,
                           99,
                           3,
                           6,
                           99,
                           99,
                           4,
                           3,
                           3,
                           3,
                           3,
                           3,
                           3,
                           99,
                           3,
                           3,
                           3,
                           99,
                           4,
                           6,
                           3,
                           3,
                           99,
                           3,
                           3,
                           3,
                           8,
                           3,
                           99,
                           99,
                           99,
                           6,
                           4,
                           3)
#done this would be the major grouping of them
length(which(sig_pathway_groupings==3))
length(which(sig_pathway_groupings==4))
length(which(sig_pathway_groupings==6))


#add type data to the sig pathways dataframe 
sig_pathways_grp <- rbind.data.frame(sig_pathways, sig_pathway_groupings)

###function that gets the mean differences between samples from IC and NIC patients
get_mean_diff <- function(x){
  
  #quit if rownames for metadata don't match input
  #stopifnot(rownames(metadata_fix) == rownames(x))
  
  #get infection means
  infection_means <- colMeans(x[which(metadata_fix$Infection_in_6=="Y"),])
  
  #non infection means
  noninfection_means <- colMeans(x[which(metadata_fix$Infection_in_6=="No"),])
  
  #Mean Dif
  Mean_Dif <- noninfection_means - infection_means
  
  #bind to new dataframe
  Pathways <- gsub(".*:","",colnames(x))
  Pathways <- gsub("superpathway of", "", Pathways)
  ret_frame <- rbind.data.frame(x, Mean_Dif, Pathways, stringsAsFactors = F)
  rownames(ret_frame)[nrow(ret_frame)] <- "Pathways"
  rownames(ret_frame)[nrow(ret_frame)-1] <- "Mean_Diff"
  ret_frame <- as.data.frame(t(ret_frame))
  ret_frame$Mean_Diff <- as.numeric(as.character(ret_frame$Mean_Diff))
  return(ret_frame)
}

all_mean_diff <- get_mean_diff(sig_pathways_grp)

#make sure metadata and sig_pathways was in same order so function works correctly.
rownames(sig_pathways_grp) == rownames(metadata_fix)

#function that converts type codes to text
to_writting <- function(x){
  for(i in 1: nrow(x) ){
    if(x$'45'[i]=='3'){
      print(rownames(x)[i])
      x$Type[i] <- "Biosynthetic"
      x$Color[i] <- "Orange"
    }else if(x$'45'[i]=='4'){
      x$Type[i]<-"Degradation"
      x$Color[i] <- "green"
    }else if(x$'45'[i]=='99'){
      x$Type[i]<-"Remove"
      x$Color[i] <- "Black"
    }else{
      x$Type[i]<-"Other"
      x$Color[i] <- "Pink"}
  }
  return(x)
}


all_mean_diff_test <- to_writting(all_mean_diff)
#remove ones that need to be removed (ones that were marked as 99 earlier)
all_mean_filt <- all_mean_diff_test[which(all_mean_diff_test$Type != "Remove"),]



#generate plot to see mean differences between two groups 
#(not included in manuscript but is similar to plot 4 without the breakdown by genera)
library(ggplot2)
library(gridGraphics)
#make base of plot
test3 <-ggplot(all_mean_filt, aes(x=Pathways, y=Mean_Diff, fill=Mean_Diff > 0))
#add geoms
test3 <- test3 + geom_rect(data=all_mean_filt, aes(x=NULL, y=NULL, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=Color), alpha=1) + 
  geom_col(position="dodge") + 
  scale_fill_manual(values=c("Red","green","Orange","Pink","blue")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=18 ), axis.text.y=element_text(size=10), axis.title.x = element_text(size=18)) +
  xlab("") + scale_y_continuous(labels = scales::percent) + ylab("Mean Difference in Relative Abundance") +
  guides(fill=F) +
  facet_grid(rows=vars(all_mean_filt$Type), scales="free_y", space="free_y", drop=T) + 
  coord_flip()
test3
#color facet labels
test2 <- ggplot_gtable(ggplot_build(test3))
stripr <- which(grepl('strip-r', test2$layout$name))
fills <- c("Orange","green","Pink")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', test2$grobs[[i]]$grobs[[1]]$childrenOrder))
  test2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.newpage()
grid.draw(test2)


#write function that will be used in future contribution scripts
write.table(all_mean_filt, file="~/projects/MALL-MGS-ALL/Sig_Pathways.tsv", col.names=T, row.names=T, quote=F, sep="\t")




















###biosynth pathways
sig_biosynth <- sig_pathways_grp[-nrow(sig_pathways_grp),which(sig_pathways_grp["45",]==3)]

####degradation
sig_degrad <- sig_pathways_grp[-nrow(sig_pathways_grp),which(sig_pathways_grp["45",]==4)]
####metabolism
sig_energy <- sig_pathways_grp[-nrow(sig_pathways_grp),which(sig_pathways_grp["45",]==6)]
#### Others
sig_others <- sig_pathways_grp[-nrow(sig_pathways_grp),which(sig_pathways_grp["45",] != 3 & sig_pathways_grp["45",] != 4 & sig_pathways_grp["45",] != 99)]

#Now we have them all split
#lets write a function that gets the mean difference between the two groups
#takes in a data.frame

sig_biosynth_means <- get_mean_diff(sig_biosynth)
sig_degrad_means <- get_mean_diff(sig_degrad)
sig_energy_means <- get_mean_diff(sig_energy)
sig_others_means <- get_mean_diff(sig_others)

#plot the mean differences for each

#Biosynth plot remember that they are in percentages
library(ggplot2)

#Biosynth plot
biosynth_plot <- ggplot(sig_biosynth_means, aes(Pathways, Mean_Diff)) + geom_bar(stat="identity", aes(fill=Mean_Diff > 0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=18, color="black" ), axis.text.y=element_text(size=10, color="black"), axis.title.x = element_text(size=18)) +
  xlab("") + scale_y_continuous(labels = scales::percent) + coord_flip() + guides(fill=F) + ggtitle("Biosynthetic Pathways")
biosynth_plot

#degrad plot
degrad_plot <- ggplot(sig_degrad_means, aes(Pathways, Mean_Diff)) + geom_bar(stat="identity", aes(fill=Mean_Diff > 0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=18, color="black" ), axis.text.y=element_text(size=10, color="black"), axis.title.x = element_text(size=18)) +
  xlab("") + scale_y_continuous(labels = scales::percent) + coord_flip() + ggtitle("Degradation Pathways") + 
  scale_fill_discrete(name="Difference in RA", breaks=c("TRUE","FALSE"),labels=c("Higher in NIC", "Higher in IC"))
degrad_plot

#energy metabolism plot
#going to combine this and other together.
#don't need this table anymore
energy_plot <- ggplot(sig_energy_means, aes(Pathways, Mean_Diff)) + geom_bar(stat="identity", aes(fill=Mean_Diff > 0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=18, color="black" ), axis.text.y=element_text(size=10, color="black"), axis.title.x = element_text(size=18)) +
  xlab("") + scale_y_continuous(labels = scales::percent) + coord_flip() + guides(fill=F) + ggtitle("Energy Metabolism Pathways")
energy_plot

#Other plot
Other_plot <- ggplot(sig_others_means, aes(Pathways, Mean_Diff)) + geom_bar(stat="identity", aes(fill=Mean_Diff > 0)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=18, color="black" ), axis.text.y=element_text(size=10, color="black"), axis.title.x = element_text(size=18)) +
  xlab("") + scale_y_continuous(labels = scales::percent) + coord_flip() + guides(fill=F) + ggtitle("Other Pathways")
Other_plot


#just the guide
legend_test <- get_legend(degrad_plot + theme(legend.direction = "vertical", legend.justification = "center" , legend.box.just = "top"))
degrad_plot <- degrad_plot + theme(legend.position = 'none')







#we could try gene abundance
gene_table <- read.table("~/projects/MALL-MGS-ALL/humann2_final_out/humann2_genefamilies_relab_unstratified.tsv", sep="\t", header=T, row.names=1, comment.char = "", quote="")
colnames(gene_table) <- gsub("_Abundance.RPKs","",colnames(gene_table))
head(gene_table)
gene_table_flip <- as.data.frame(t(gene_table))
gene_table_fix <- gene_table_flip[rownames(metadata_fix),]
pvals_gene <- apply(gene_table_fix, 2, function(x) wilcox.test(x ~ metadata_fix$Infection_in_6, exact=F)$p.value)
qvals_gene <- p.adjust(pvals_gene, "fdr")
hits_gene <- qvals_gene < 0.05
which(hits_gene)
#thats a lot of hits :0 2458 in total... what the hell am I suppose to do with this????
sig_gene <- gene_table_fix[,which(hits_gene)]
means <- apply(sig_gene, 2, function(x) t.test(x ~ metadata_fix$Infection_in_6)$estimate)
means <- as.data.frame(means)
total_diff <- apply(means, 2, function(x) x[1]-x[2])
test <- rbind.data.frame(means, total_diff)
test_flip <- as.data.frame(t(test))
test_flip[order(3),]
