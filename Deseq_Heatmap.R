library(ggplot2)
library(DESeq2) #BiocManager::install("DESeq2") 
library(phyloseq)
library(dplyr)
library(ggplotify)

###### 
#File import 
physeq_f = readRDS("IBS.RDS")
physeq_f = tax_glom(physeq_f, "Genus")
Sample_data = read.csv("IBS_1Year.csv")
Sample_data$abp_Combined = ifelse(is.na(Sample_data$abp_Multiple) == T, NA, 
                             ifelse(Sample_data$abp_Multiple == "No", "None", "Any"))
Sample_data$Any = ifelse(Sample_data$Any == "No", "None", "Any")

taxa_identifiers = read.csv("~/Desktop/ABIS/RDP_Taxa_identifiers.csv")
row.names(taxa_identifiers) = taxa_identifiers$ASV
taxa_identifiers = taxa_identifiers[,c("ASV", "Genus")]

#################################################################
#User Input 
alpha = 0.05

#################################################################
#Abundance
#################################################################
ps.RA = transform_sample_counts(physeq_f, function(x) x / sum(x) )

qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
rownames(qpcr_data) = qpcr_data$ID
qpcr_data$ID = NULL

otu_RA = data.frame(otu_table(ps.RA))
otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
row.names(otu_RA_1) = otu_RA_1$Row.names
otu_RA_1$Row.names = NULL

for(i in 2:length(names(otu_RA_1))) {
  otu_RA_1[,i] <- round(otu_RA_1[,1] * otu_RA_1[, i])
}

otu_RA_1$copies_16s_per_gram_stool = NULL

ps.RA.qpcr = ps.RA
otu_table(ps.RA.qpcr) = otu_table(otu_RA_1, taxa_are_rows = FALSE)
ps.RA.qpcr

Actual_sign = data.frame()
current_column = "Any"
for (current_column in c("Any", "abp_Combined", "k_Combined", "r_Combined")) {
  print(current_column)
  ##################################
  #Any Vs None
  ps.current=ps.RA.qpcr
  Sample_data$Group = Sample_data[,current_column]
  rownames(Sample_data) = Sample_data$ID
  sample_data(ps.current) = sample_data(Sample_data)
  
  subset_cases = subset(Sample_data, Sample_data[,current_column] != "None" &
                          Sample_data[,current_column] != "No" )
  subset_controls = subset(Sample_data, Sample_data$Any == "None" )
  
  
  OTU_table = data.frame(otu_table(ps.current))
  OTU_table$ID = rownames(OTU_table)
  #Prevalence filter
  ASV_list = colnames(OTU_table)
  
  current_tax = ASV_list[1]
  for (current_tax in ASV_list) {
    current_subset = OTU_table[,c("ID", current_tax)]
    IBS = subset(current_subset, current_subset$ID %in% subset_cases$ID)
    IBS_prev = round(nrow(subset(IBS, IBS[,current_tax] > 0))/nrow(IBS),3)
    Controls = subset(current_subset, ! current_subset$ID %in% subset_cases$ID)
    Controls_prev = round(nrow(subset(Controls, Controls[,current_tax] > 0))/nrow(Controls),3)
    if (IBS_prev < 0.25) {
      OTU_table[,current_tax] = NULL
    }
  }
  
  OTU_table[is.na(OTU_table)] = 0
  OTU_table <- mutate_all(OTU_table, function(x) as.integer(as.character(x)))
  otu_table(ps.current) = otu_table(OTU_table, taxa_are_rows = F) #[ 141 taxa and 560 samples ]
  
  #DESEQ
  diagdds = suppressWarnings(phyloseq_to_deseq2(ps.current, ~ Group))
  diagdds = estimateSizeFactors(diagdds, type = "poscounts")
  diagdds = DESeq(diagdds, test="Wald", fitType="local", quiet = TRUE)
  res = results(diagdds, cooksCutoff = FALSE)
  print(resultsNames(diagdds))
  
  sigtab = data.frame(res)#[which(res$padj <= 0.1), ])
  sigtab$ASV = rownames(sigtab)
  sigtab= merge(sigtab, taxa_identifiers)
  
  if (nrow(sigtab)> 0) {
    sigtab[,"Genus"] = gsub("-", ".", sigtab[,"Genus"])
    Identifier_list = unique(sigtab[,"Genus"])
    significant_data = sigtab[,c("log2FoldChange", "Genus", "pvalue", "padj")]
    significant_data$Group = current_column
    Actual_sign = rbind(Actual_sign, significant_data)
  } else {
    print(paste(current_column, "No Significant"))
  }
}
Actual_sign$Genus = gsub("_", " ", Actual_sign$Genus)
Actual_sign$Genus = gsub("Escherichia/Shigella", "Escherichia/Shigella", Actual_sign$Genus)
Actual_sign_1 = Actual_sign
#Actual_sign = subset(Actual_sign, Actual_sign$Group %in% c("Any", "Multiple"))
set.seed(15)

Actual_sign$Group = gsub("_Combined", "", Actual_sign$Group)
Actual_sign$Group = gsub("Any", "AKR", Actual_sign$Group)
Actual_sign$Group  = toupper(Actual_sign$Group)
  

##########
#Pvalue Subset 
Pvalue_df = Actual_sign
Pvalue_df_1 = subset(Pvalue_df, Pvalue_df$padj <= 0.05)


Pvalue_df = subset(Actual_sign, Actual_sign$Genus %in% Pvalue_df_1$Genus)
Pvalue_df$P_type = ifelse(Pvalue_df$padj <= 0.05, "***",
                          ifelse(Pvalue_df$padj <= 0.1, "**", ""))
Pvalue_df = Pvalue_df[,c("Genus", "Group", "P_type")]
Pvalue_df = reshape(Pvalue_df, idvar = "Genus", timevar = "Group", 
                           direction = "wide")
rownames(Pvalue_df) = Pvalue_df$Genus; Pvalue_df$Genus = NULL
colnames(Pvalue_df) = gsub("P_type.", "", colnames(Pvalue_df))
Pvalue_df[is.na(Pvalue_df)] = ""
Pvalue_df = Pvalue_df[,c("AKR", "ABP", "K", "R")]
Pvalue_df = Pvalue_df[order(rownames(Pvalue_df)),]

Actual_sign_long = subset(Actual_sign, Actual_sign$Genus %in% rownames(Pvalue_df))
Actual_sign_long= Actual_sign_long[,c("Genus", "Group", "log2FoldChange")]
Actual_sign_long = reshape(Actual_sign_long, idvar = "Genus", timevar = "Group", 
                           direction = "wide")
rownames(Actual_sign_long) = Actual_sign_long$Genus; Actual_sign_long$Genus = NULL
colnames(Actual_sign_long) = gsub("log2FoldChange.", "", colnames(Actual_sign_long))
Actual_sign_long = Actual_sign_long[order(rownames(Actual_sign_long)),]


cols = c(colorRampPalette(c("darkblue", "lightskyblue"))(100),
         colorRampPalette(c("brown4","chocolate1"))(100))

breaks = c(seq(min(Actual_sign_long,na.rm = T),-0.0001,length.out = 100), 0,
           seq(max(Actual_sign_long, na.rm = T),0.0001,length.out = 100))

#Actual_sign_long[is.na(Actual_sign_long)] = 0

Actual_sign_long = Actual_sign_long[order(rownames(Actual_sign_long)),]

heatmap = pheatmap(Actual_sign_long, 
                   cluster_cols = F,
                   cluster_rows = T,
                   color = cols, breaks = breaks,
                   display_numbers = Pvalue_df,
                   number_color = "black",
                   fontsize_number = 12,
                   treeheight_row = 25,
                   treeheight_col = 15,
                   na_col = "grey65", border_color = "black")


jpeg("./Images/Deseq_Groups_heatmap.jpeg", 
     res = 500, height = 2000, width = 4000)
as.ggplot(heatmap)
dev.off()
