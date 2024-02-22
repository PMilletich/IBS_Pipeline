#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(ggplot2)
library(metagMisc)

#################################
###File Import 
setwd("~/Desktop/IBS/")

set.seed(1)

physeq_f = readRDS("./IBS.RDS")
physeq_f = tax_glom(physeq_f, "Phylum")
sample = data.frame(sample_data(physeq_f))

#################################################################
#Relative Abundance
#################################################################
ps.RA = transform_sample_counts(physeq_f, function(x) x / sum(x) )

#################################################################
#Create long PS 
#################################################################
OG_PS = ps.RA
OG_taxa = tax_table(OG_PS)

ibs_col = c("Control", "cdsrABP",  "cdABP", "srABP")

long_sample = data.frame()
long_otu = data.frame()

for(current_col in ibs_col) {
  if (current_col == "Control") {
    current_subset = subset(sample, sample$IBS == "Control")
  } else if (current_col == "cdsrABP") {
    current_subset = subset(sample, sample$IBS != "Control" )
  } else {
    current_subset = subset(sample, sample$IBS == current_col)
  } 
  
  IBS_df = unique(current_subset)
  rownames(IBS_df) = IBS_df$ID
  IBS_df$Group = current_col
  print(table(IBS_df$Group))
  
  current_PS = OG_PS
  sample_data(current_PS) = sample_data(IBS_df)
  #Sample
  current_sample_data = data.frame(sample_data(current_PS))
  rownames(current_sample_data) = paste(current_sample_data$ID, current_col, sep = "_")
  long_sample = rbind(long_sample, current_sample_data)
  #OTU
  current_otu_data = data.frame(otu_table(current_PS))
  rownames(current_otu_data) = paste(rownames(current_otu_data), current_col, sep = "_")
  long_otu = rbind(long_otu, current_otu_data)
}

table(long_sample$Group )

long_ps = phyloseq(otu_table(long_otu, taxa_are_rows = F), 
                   tax_table(as.matrix(OG_taxa)), 
                   sample_data(long_sample))


#########################################################
#Create Master Table
current_ps = long_ps

current_ps = phyloseq_rename_with_tax(current_ps, "Phylum")
count_table = data.frame(otu_table(current_ps))
taxa_list = colnames(count_table)
count_table$ID = rownames(count_table)

long_sample$ID = rownames(long_sample)
Total_data = merge(count_table, long_sample)

###############################################################
#Relative abundance of Phylum differences between IBS groups 
###############################################################

Abundance_df = data.frame()
Prevalence_df = data.frame()
current_taxa = "Firmicutes"
for (current_taxa in taxa_list) {
  current_data = Total_data
  Prevalence = round(nrow(subset(current_data, current_data[,current_taxa] > 0))/ nrow(current_data),3)*100
  Prevalence_row = data.frame("Prev" = Prevalence, 
                              "Taxa" = current_taxa)
  Prevalence_df = rbind(Prevalence_df, Prevalence_row)
  
  for (current_group in unique(current_data$Group)) {
    current_data_subset = subset(current_data, current_data$Group == current_group)
    Abund_row = data.frame("Taxa" = current_taxa, 
                           "Group" = current_group,
                           "Abundance" = current_data_subset[,current_taxa])
    Abundance_df = rbind(Abundance_df, Abund_row)
  }
}

###############################################################
#Firmicutes/Bacteroides Ratio 
###############################################################
current_data = Total_data

for (current_group in unique(current_data$Group)) {
  current_data_subset = subset(current_data, current_data$Group == current_group)
  current_data_subset$FB_Ratio = current_data_subset$Firmicutes/current_data_subset$Bacteroidetes
  xx = current_data_subset[,c("ID", "Firmicutes", "Bacteroidetes", "FB_Ratio")]
  
  Abund_row = data.frame("Taxa" = "FB_Ratio", 
                         "Group" = current_group,
                         "Abundance" = current_data_subset$FB_Ratio)
  Abundance_df = rbind(Abundance_df, Abund_row)
}

###############################################################
#Prevalence thresholding
###############################################################
Prevalence_df = Prevalence_df[order(Prevalence_df$Prev, decreasing = T),]
Prevalence_df = subset(Prevalence_df, Prevalence_df$Prev >= 50)
Abundance_df_subset = subset(Abundance_df, Abundance_df$Taxa %in% Prevalence_df$Taxa)

###############################################################
#Reorder X axis and create list of comparisons for stat_compare_means
###############################################################
Abundance_df_subset$Group = factor(Abundance_df_subset$Group, 
                                   levels = c("Control","cdsrABP",
                                   "srABP", 
                                   "cdABP"))
comparisons_list= list(c("cdsrABP","Control"),
                       c("srABP","Control"), 
                       c("cdABP","Control"))

###############################################################
#Create and save image 
###############################################################
Abundance_df_subset$RelativeAbundance = Abundance_df_subset$Abundance*100
All_Phyla = ggplot(Abundance_df_subset, aes(x = Group,
                                            y = RelativeAbundance, 
                                            fill = Group,
                                            color = Group)) + 
  geom_boxplot() + theme_bw() + theme(legend.position = "top") + 
  facet_wrap(~Taxa, scales = "free", ncol = 3) + 
  theme(axis.text.x = element_text(size = 8), 
        legend.position = "none",
        axis.title.x = element_blank()) + 
  stat_compare_means(comparisons = comparisons_list, 
                     step.increase = 0.075)+ 
  scale_fill_manual(name = "IBS Group",
                    breaks = c("cdsrABP", "srABP", 
                               "cdABP", "Control"), 
                    values = c("#65ACEA", "#F9D057",
                               "#D05783", "grey25"))+ 
  scale_color_manual(name = "IBS Groups", 
                     breaks = c("cdsrABP", "srABP",
                                "cdABP", "Control"), 
                     values = c("darkblue", "goldenrod",
                                "deeppink4", "grey5")); All_Phyla

jpeg("./Images/IBS_4/Phylum_3.jpeg", 
     res = 500, height = 4000, width = 4000)
All_Phyla
dev.off()
