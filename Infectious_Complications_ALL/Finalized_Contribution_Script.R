#Cleaned up Script for Contribution breakdown figure.


#read in metadata
metadata <- read.table("~/projects/MALL-deblur/metadata/After_updates_to_redacap/metdata_updated_from_July_17.csv",
                       sep="\t", header=T, row.names = 1, stringsAsFactors = F, comment.char="")

#read in table containing the significantly different pathways

sig_pathways <- read.table("~/projects/MALL-MGS-ALL/Sig_Pathways.tsv", sep="\t", header=T, row.names = 1)

#load in stratified data

stratified_data <- read.table("~/projects/MALL-MGS-ALL/humann2_final_out/humann2_pathabundance_relab_stratified.tsv",
                              sep="\t", header=T, row.names=1, comment.char="", quote="")

#okay now need to fix up metadata.
#fix rownames
rownames(metadata) <- gsub("-",".",rownames(metadata))
#only keep metadata on samples that are in 6 months (therefore contained in sig_pathways table)
metadata_fix <- metadata[colnames(sig_pathways[1:44]),]


#break strat data into only sig_strat pathways
#fix col names
colnames(stratified_data) <- gsub("_Abundance","",colnames(stratified_data))

#only keep strat samples in 6 months
stratified_data_fix <- stratified_data[,rownames(metadata_fix)]

#okay now we can only keep those that are sig different.

#make empty data frame
sig_strat_fix_data <- data.frame()

#get pathway codes for each significant different pathway and then pull them out of the stratified table.
for(i in rownames(sig_pathways)){
  #pull out code
  Pathway <- gsub(":.*",":",i)
  #check if the code matchs first thing in the string name
  Pathway <- paste0("^",Pathway)
  #add to data frame
  sig_strat_fix_data <- rbind(sig_strat_fix_data, stratified_data_fix[grep(Pathway, rownames(stratified_data_fix)),])
}

#okay awesome now we need to get it in a form that we can make graphs with in in ggplot.

#add pathway name coloumn to each row
#remove pathway code names
pathway_names <- gsub(".*:","",rownames(sig_strat_fix_data))
pathway_names
#remove bacteria names
pathway_names <- gsub("\\|.*","",pathway_names)
pathway_names
#remove superpathway due to it taking a lot of space
pathway_names <- gsub("superpathway of", "", pathway_names)
pathway_names
#add to sig_start data
sig_strat_fix_data$Pathway <- pathway_names

#okay we need the species name and the genus name.
#pull out bacteria names
bacteria <- gsub(".*\\|","",rownames(sig_strat_fix_data))
bacteria
#get spceices names
bacteria <- gsub(".*s__","",bacteria)
bacteria

#set them
sig_strat_fix_data$Species <- bacteria

#get genus
genus <- gsub("_.*","", bacteria)
genus

#set genus name
sig_strat_fix_data$genus <- genus
#set the type of pathway each one is

sig_strat_fix_data$Type <- sig_pathways[match(sig_strat_fix_data$Pathway, sig_pathways$Pathways), "Type"]



#okay now we want to melt the data


library(reshape2)

sig_strat_melt <- melt(sig_strat_fix_data, id.vars = c("Species", "genus", "Type", "Pathway"))

#make an infectious complciation variable
sig_strat_melt$Infection <- metadata_fix[sig_strat_melt$variable, "Infection_in_6"]
#change Y to Yes
sig_strat_melt$Infection <- ifelse(sig_strat_melt$Infection=="Y", "Yes", "No")

#create genus collapse column
sig_strat_melt$gen_col <- as.character(sig_strat_melt$genus)
#okay now we need make another genus column thats collapse if the total RA is below 0.001
collapse_genus <- function(x){
  for (bac in unique(x$genus)){
    total = 0
    for (i in 1:nrow(x)){
      if(x$genus[i] == bac){
        total = total + x$value[i]
      }
    
    }
    print(total)
    if(total <= 0.0001){
      print("replacing genus")
      for(i in 1:nrow(x)){
        if(x$genus[i] == bac){
          x$gen_col[i] <- "zOther"
        }
      }
    }
  }
  return(x)
}


remove_zero_contributors <- function(x){
  for (path in unique(x$Pathway)){
    total = 0
    for (i in 1:nrow(x)){
      if(x$Pathway[i] == path){
        total = total + x$value[i]
      }
    }
    if(total <= 0.00001){
      print("removing pathway")
      print(path)
      for(i in 1:nrow(x)){
        if(x$Pathway[i] == path){
          x$keep[i] <- "remove"
        }else
          x$keep[i] <- "keep"
      }
    }
  }
  return(x)
  
}
#collapse genus
sig_strat_melt_col <- collapse_genus(sig_strat_melt)
#remove pathways that sum to zero due to no bugs continuing the whole rxn
sig_strat_melt_col <- remove_zero_contributors(sig_strat_melt_col)

sig_strat_melt_col<- subset(sig_strat_melt_col, subset=keep=="keep")

#okay the data is all prepared we can now go ahead and plot it :^]
library(ggplot2)
#contributation plot
#get color pallet.
diff_col <- c("#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
              "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
              "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
              "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
              "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
              "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
              "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
              "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
              "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
              "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
              "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
              "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
              "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

metled_plot <- ggplot(sig_strat_melt_col, aes(Infection, value, fill=gen_col)) + geom_bar(stat="summary", ) +
  guides(fill=guide_legend(nrow=51, title="Genus")) + theme(legend.key.size = unit(.75, "line"), strip.text.y = element_text(angle=180, hjust=1),
                                                           panel.grid.major = element_blank(),
                                                           panel.grid.minor = element_blank(),
                                                           panel.border = element_blank(),
                                                           panel.background = element_blank(),
                                                           strip.background = element_blank(),
                                                           panel.spacing = unit(.05,"lines"),
                                                           strip.placement = "outside") +
  ylab("Mean Relative Abundance") +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  facet_grid(rows=vars(sig_strat_melt_col$Pathway), switch="y") +
  scale_fill_manual(values=diff_col) +
  coord_flip()

metled_plot

#we can make a similar plot but we want back to back bars

test_plot <- ggplot(sig_strat_melt_col, aes(x=Pathway, value, fill=gen_col)) +
  geom_bar(data=subset(sig_strat_melt_col, Infection=="Yes"), aes(y=value), stat="summary") +
  geom_bar(data=subset(sig_strat_melt_col, Infection=="No"), aes(y=-value), stat="summary") +
  geom_hline(yintercept = 0) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=10, colour = "black") ) +
  scale_fill_manual(values=diff_col) +
  guides(fill=guide_legend(title="Genus")) +
  ylab("Mean Relative Abundance") +
  facet_grid(Type ~ ., space="free", scales="free_y") +
  coord_flip() +
  scale_y_continuous(labels=scales::percent)
test_plot

#color the facet labels
final_plot <- ggplot_gtable(ggplot_build(test_plot))
stripr <- which(grepl('strip-r', final_plot$layout$name))
fills <- c("#e2df0f","#eac4ed","Orange")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', final_plot$grobs[[i]]$grobs[[1]]$childrenOrder))
  final_plot$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
library(gridGraphics)
grid.newpage()
grid.draw(final_plot)

