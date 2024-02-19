###############################################################
#Load Files and packages
library(ggplot2) #https://ggplot2.tidyverse.org/
library(pheatmap) #https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap
#library(ggplotify)

setwd("~/Desktop/IBS/")
OG_PS = readRDS("IBS.RDS")
OG_PS = tax_glom(OG_PS, "Genus")
OG_taxa = data.frame(tax_table(OG_PS))
map = read.csv("IBS_1Year_IBS.csv")

##############
#Import significant genera from pime, lefse and aldex
PIME = read.csv("./Images/IBS_4/PIME_All_4_1.csv")
Aldex = read.csv("./Images/IBS_4/Aldex2.csv")
Lefse = read.csv("./Images/IBS_4/Lefse.csv")

length(unique(PIME$Genus)) #44
length(unique(Aldex$feature)) #16
length(unique(Lefse$feature)) #14

#Find genera overlapping between P&A P&L
PIME_Lefse = intersect(PIME$Genus, Lefse$feature)
PIME_Aldex = intersect(PIME$Genus, Aldex$feature)
MDA_Count = unique(c(PIME_Aldex, PIME_Lefse))

#Subset if TAXA were found to be significant 
OG_taxa = subset(OG_taxa, OG_taxa$Genus %in% MDA_Count)

#################################################################
#Abundance 
#################################################################
ps.RA = transform_sample_counts(OG_PS, function(x) x / sum(x) )

qpcr_data = map[,c("ID", "copies_16s_per_gram_stool")]
rownames(qpcr_data) = qpcr_data$ID
otu_RA = data.frame(otu_table(ps.RA))
otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
row.names(otu_RA_1) = otu_RA_1$Row.names
otu_RA_1$Row.names = NULL
otu_RA_1$ID = NULL
otu_RA_2 = otu_RA_1

for(i in 2:length(names(otu_RA_2))) {
  otu_RA_2[,i] <- round(otu_RA_2[,1] * otu_RA_2[, i])
}

ps.qpcr = ps.RA
otu_table(ps.qpcr) = otu_table(otu_RA_2, taxa_are_rows = F)

#Create total dataframe of genera and metadata
current_ps = ps.qpcr
current_ps = phyloseq_rename_with_tax(current_ps, "Genus")
otu_df = data.frame(otu_table(current_ps))
otu_df$ID = rownames(otu_df)
Total_data = merge(otu_df, map)

###############
#For each genera, find the difference in control and abp/ibs mean
Group_data = data.frame()
logchange_data = data.frame()

for (current_genera in unique(MDA_Count)) {
  print(current_genera)
  Current_Sign = data.frame()
  current_genera = gsub("/", ".", current_genera)
  for (current_Group in c("cdsrABP", "srABP", "cdABP")) {
    if (current_Group == "cdsrABP") {
      IBS_df = subset(Total_data, Total_data$IBS != "Control")
    } else {
      IBS_df = subset(Total_data, Total_data$IBS == current_Group)
    }
    
    Control = subset(Total_data, Total_data$IBS == "Control")
    print(paste(current_Group, nrow(Control), nrow(IBS_df)))
    
    LogChange = log2(mean(Control[,current_genera], na.rm = T)/mean(IBS_df[,current_genera]))

    current_row = data.frame("Taxa" = current_genera,
                             "Group" = current_Group,
                             "LogChange" = LogChange)
    Current_Sign = rbind(Current_Sign, current_row)
    
    #Determine label if its in PIME, LEFSE and/or aldex
    Group_Index = ""

    PIME_subset = subset(PIME, PIME$Genus == gsub("\\.", "/", current_genera) & 
                           PIME$Group == current_Group)
    
    if (nrow(PIME_subset) ==1) {
      Group_Index =paste(Group_Index, "P", sep = "")
    }
    
    Lefse_subset = subset(Lefse, Lefse$feature ==  gsub("\\.", "/", current_genera) & 
                            Lefse$Group == current_Group)
    
    if (nrow(Lefse_subset) ==1) {
      Group_Index =paste(Group_Index, "L", sep = "")
    }
    
    Aldex_subset = subset(Aldex, Aldex$feature ==  gsub("\\.", "/", current_genera) & 
                            Aldex$Group == current_Group)
    if (nrow(Aldex_subset) ==1) {
      Group_Index =paste(Group_Index, "A", sep = "")
    }
    
    Group_data= rbind(Group_data, 
                      data.frame("Taxa" = current_genera, 
                                 "Group" = current_Group,
                      "Index" = Group_Index))
  }
  logchange_data = rbind(logchange_data, Current_Sign)
}

logchange_data$Taxa = gsub("Escherichia.Shigella", 
                           "Escherichia/Shigella", logchange_data$Taxa )
Group_data$Taxa = gsub("Escherichia.Shigella", 
                           "Escherichia/Shigella", Group_data$Taxa )
########################################################
#convert dataframe to matrix for pheatmap
log_data = logchange_data[,c("Taxa", "Group", "LogChange")]
log_data_1= reshape(log_data, idvar = "Taxa", timevar = "Group", direction = "wide")
rownames(log_data_1) = log_data_1$Taxa; log_data_1$Taxa = NULL
colnames(log_data_1) = gsub("LogChange.", "", colnames(log_data_1))


Index_data = reshape(Group_data, idvar = "Taxa",
                     timevar = "Group", direction = "wide")
colnames(Index_data) = gsub("Index.", "", colnames(Index_data))
rownames(Index_data) = Index_data$Taxa; Index_data$Taxa = NULL

#ensure columns are in the same order
log_data_1 = log_data_1[c("cdsrABP", "srABP", "cdABP")]
Index_data = Index_data[c("cdsrABP", "srABP", "cdABP")]

#Create annotations for phylum
tax_df = data.frame(tax_table(current_ps))

Phylum = tax_df[,c("Phylum", "Genus")]
rownames(Phylum) = Phylum$Genus
Phylum$Genus = NULL
ann_colors = list(Phylum = c("Actinobacteria"="darkgreen", 
                             "Firmicutes"="goldenrod1",
                             "Proteobacteria"="orchid3",
                             "Bacteroidetes"="purple4"))

#Create manual color scale
breaks = c(seq(min(log_data_1),-0.0001,length.out = 100), 0,
           seq(max(log_data_1),0.0001,length.out = 100))

cols = c(colorRampPalette(c("brown4","sienna1"))(100),
         colorRampPalette(c("cornflowerblue", "lightblue1"))(100))


heatmap = pheatmap(log_data_1, 
                   cluster_cols = F,
                   cluster_rows = T,
                   annotation_row = Phylum,
                   annotation_colors = ann_colors,
                   color = cols, breaks = breaks,
                   display_numbers = Index_data,
                   number_color = "black",
                   fontsize_number = 10,
                   treeheight_row = 25,
                   treeheight_col = 15,
                   na_col = "grey15", border_color = "black")

jpeg("./Images/IBS_4/PIME_heatmap_all.jpeg", res = 600, 
     height = 3000, width = 4000)
heatmap
dev.off()
