#Create a long dataframe of the IBS data 
#####################
#Load Files 
#####################
library(haven)
library(ggplot2) #https://ggplot2.tidyverse.org/
library(ggpubr) 
library(pheatmap)
library(ggplotify)


setwd("~/Desktop/IBS/")
#Sample Data and Phyloseq from Previous Subsetting Step 
OG_PS = readRDS("IBS.RDS")
OG_PS = tax_glom(OG_PS, "Genus")
OG_taxa = data.frame(tax_table(OG_PS))
map = read.csv("IBS_1Year.csv")

Genera = read.csv("./CSV_Files/PIME_All.csv")
MDA_Count = data.frame(table(Genera$Genus))
length(MDA_Count$Var1)
MDA_Count = subset(MDA_Count, MDA_Count$Freq > 1)
OG_taxa = subset(OG_taxa, OG_taxa$Genus %in% MDA_Count$Var1)

#################################################################
#Relative Abundance 
#################################################################
ps.RA = transform_sample_counts(OG_PS, function(x) x / sum(x) )

#################################################################
#qPCR 
#################################################################
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

#############
#Long data and new Phyloseq
OG_OTU = data.frame(otu_table(ps.qpcr))

long_sample = data.frame()
long_otu = data.frame()

ibs_col = c("abp_Multiple","k_Combined","r_Combined","Any", "None")
current_col = ibs_col[1]
for(current_col in ibs_col) {
  if (current_col != "None") {
    current_subset = subset(map, map[,current_col] != "No" &  
                              map[,current_col] != "None")
  } else {
    current_subset = subset(map, map$Any == "No")
  }
  print(paste(current_col, nrow(current_subset)))
  rownames(current_subset) = current_subset$ID
  current_PS = OG_PS
  sample_data(current_PS) = sample_data(current_subset)
  #Sample
  current_sample_data = data.frame(sample_data(current_PS))
  rownames(current_sample_data) = paste(current_sample_data$ID, current_col, sep = "_")
  current_sample_data$Group = current_col
  long_sample = rbind(long_sample, current_sample_data)
  #OTU
  current_otu_data = data.frame(otu_table(current_PS))
  rownames(current_otu_data) = paste(rownames(current_otu_data), current_col, sep = "_")
  long_otu = rbind(long_otu, current_otu_data)
}
long_sample$Group = ifelse(long_sample$Group == "abp_Multiple", "abp_Combined", long_sample$Group)
long_sample$Group = factor(long_sample$Group, 
                           levels = c("None","Any", 
                                      "abp_Combined","k_Combined","r_Combined"))
table(long_sample$Group)
long_ps = phyloseq(otu_table(long_otu, taxa_are_rows = F), 
                  tax_table(as.matrix(OG_taxa)), 
                  sample_data(long_sample))

#############
#Total Data 
long_sample$ID_Groups = rownames(long_sample)

OG_taxa$ASV = rownames(OG_taxa)
OG_taxa = OG_taxa[,c("ASV", "Genus")]

long_otu_t = data.frame(t(long_otu))
long_otu_t$ASV = rownames(long_otu_t)
long_otu_t = subset(long_otu_t, long_otu_t$ASV %in% OG_taxa$ASV)
long_otu_t= merge(long_otu_t, OG_taxa)
rownames(long_otu_t) = long_otu_t$Genus
long_otu_t$ASV = NULL; long_otu_t$Genus = NULL
long_otu = data.frame(t(long_otu_t))
rownames(long_otu) = sub('.', '', rownames(long_otu))
long_otu$ID_Groups = rownames(long_otu)

Total_data = merge(long_otu, long_sample)

Sign_data = data.frame()
Prev_data = data.frame()
for (Genera_4 in unique(MDA_Count$Var1)) {
  Current_Sign = data.frame()
  Genera_4 = gsub("/", ".", Genera_4)
  for (current_Group in c("Any", "abp_Combined", "k_Combined", "r_Combined")) {
    IBS = subset(Total_data, Total_data$Group == current_Group)
    Control = subset(Total_data, Total_data$Group == "None")
    Wilcox = wilcox.test(IBS[,Genera_4], Control[,Genera_4])
    P = Wilcox$p.value
    
    Control_Prev = round(nrow(subset(Control, Control[,Genera_4] > 0))/nrow(Control),3)*100
    IBS_Prev = round(nrow(subset(IBS, IBS[,Genera_4] > 0))/nrow(IBS),3)*100
    Prev_Diff = Control_Prev - IBS_Prev
    
    P_row = data.frame("Taxa" = Genera_4,
                       "Group" = c(current_Group, "Control"), 
                       "PrevDiff" = c(IBS_Prev,Control_Prev))
    Prev_data = rbind(Prev_data, P_row)
    LogChange = log2(mean(Control[,Genera_4], na.rm = T)/mean(IBS[,Genera_4]))

    current_row = data.frame("Taxa" = Genera_4,
                             "Group" = current_Group,
                             "PValue" = P, 
                             "LogChange" = LogChange, 
                             "PrevDiff" = Prev_Diff)
    Current_Sign = rbind(Current_Sign, current_row)
  }
  #Current_Sign$padj = p.adjust(Current_Sign$PValue)
  Sign_data = rbind(Sign_data, Current_Sign)
}

Sign_data$Taxa = gsub("Escherichia.Shigella", "Escherichia/Shigella",Sign_data$Taxa)
Sign_data$Taxa = gsub("Clostridium_sensu_stricto", "Clostridium sensu stricto",Sign_data$Taxa)
Sign_data$Group = gsub("_Combined", "", Sign_data$Group)
Sign_data$Group = gsub("Any", "AKR", Sign_data$Group)
Sign_data$Group = toupper(Sign_data$Group)


Prev_data$Group = gsub("_Combined", "", Prev_data$Group)
Prev_data$Group = gsub("Any", "AKR", Prev_data$Group)
Prev_data$Group = toupper(Prev_data$Group)
Prev_data$Group = gsub("CONTROL", "Control", Prev_data$Group)

Prev_data$Taxa = gsub("Escherichia.Shigella", "Escherichia/Shigella",Prev_data$Taxa)
Prev_data$Taxa = gsub("Clostridium_sensu_stricto", "Clostridium sensu stricto",Prev_data$Taxa)

#Prev_data$Group = toupper(Prev_data$Group)
#######
#######
#P Matrix
P_data = Sign_data[,c("Taxa", "Group", "PValue")]
P_data$padj = p.adjust(P_data$PValue)
P_data$P_type = ifelse(P_data$padj <= 0.1, "***",
                       ifelse(P_data$PValue <= 0.05, "*", ""))
P_data = subset(P_data, P_data$PValue <= 0.05)
P_data = P_data[,c("Taxa", "Group", "P_type")]

P_data= reshape(P_data, idvar = "Taxa", timevar = "Group", direction = "wide")
rownames(P_data) = P_data$Taxa; P_data$Taxa = NULL
colnames(P_data) = gsub("P_type.", "", colnames(P_data))
P_data[is.na(P_data)] = ""

log_data = Sign_data[,c("Taxa", "Group", "LogChange")]
log_data = subset(log_data, log_data$Taxa %in% rownames(P_data))
log_data_1= reshape(log_data, idvar = "Taxa", timevar = "Group", direction = "wide")
rownames(log_data_1) = log_data_1$Taxa; log_data_1$Taxa = NULL
colnames(log_data_1) = gsub("LogChange.", "", colnames(log_data_1))

cols = c(colorRampPalette(c("darkblue", "lightskyblue"))(100),
         colorRampPalette(c("brown4","chocolate1"))(100))

breaks = c(seq(min(log_data_1),-0.0001,length.out = 100), 0,
           seq(max(log_data_1),0.0001,length.out = 100))

P_data = P_data[c("AKR", "ABP", "K", "R")]
log_data_1 = log_data_1[c("AKR", "ABP", "K", "R")]

heatmap = pheatmap(log_data_1, 
                   cluster_cols = F,
                   cluster_rows = T,
                   color = cols, breaks = breaks,
                   display_numbers = P_data,
                   number_color = "black",
                   fontsize_number = 12,
                   treeheight_row = 25,
                   treeheight_col = 15,
                   #main = "log2(mean(Controls)/mean(IBS))",
                   na_col = "grey15", border_color = "black")

jpeg("./Images/PIME_heatmap.jpeg", res = 600, height = 2000, width = 4000)
heatmap
dev.off()
