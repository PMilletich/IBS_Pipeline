#Create a long dataframe of the IBS data 
#####################
#Load Files 
#####################
library(ggplot2) #https://ggplot2.tidyverse.org/
library(ggpubr) 
library(phyloseq)

#Sample Data and Phyloseq from Previous Subsetting Step 
OG_PS = readRDS("IBS.RDS")
OG_taxa = data.frame(tax_table(OG_PS))
map = read.csv("IBS_1Year.csv")

#################################################################
#Relative Abundance 
#################################################################
ps.RA = transform_sample_counts(OG_PS, function(x) x / sum(x) )

#################################################################
#qPCR; Absolute Abundance
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
long_sample$Group = gsub("_Combined", "", long_sample$Group)
long_sample$Group = gsub("None", "Controls", long_sample$Group)
long_sample$Group = gsub("Any", "AKR", long_sample$Group)


long_sample$Group = factor(long_sample$Group, 
                           levels = c("Controls","AKR", 
                                      "abp","k","r"))
table(long_sample$Group)


long_ps = phyloseq(otu_table(long_otu, taxa_are_rows = F), 
                  tax_table(as.matrix(OG_taxa)), 
                  sample_data(long_sample))

my_comparisons <- list( c("Controls", "AKR"), c("Controls", "abp"), 
                        c("Controls", "k"), c("Controls", "r"))

Physeq_rare = rarefy_even_depth(long_ps, sample.size = 25000,
                                rngseed = 13, replace = F, 
                                trimOTUs = F, verbose = T)

Alpha_rare = plot_richness(Physeq_rare, x="Group", measures=c("Observed", "Shannon"),
                      color = "Group") + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox") +  
  geom_jitter(alpha = 0.5, width = 0.4) +
  scale_color_manual(name = "IBS Groups", 
                    breaks = c("abp", "AKR", 
                               "k", "r",
                               "Controls"),
                    values = c("purple", "coral", 
                               "forestgreen", "dodgerblue", 
                               "grey")) + 
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        axis.title.x = element_blank(), 
        legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Alpha_rare$layers = Alpha_rare$layers[-1]
Alpha_rare
                                
jpeg("./Images/Alpha_rare_All_IBS.jpeg", res = 600, 
     width = 6000, height = 3000)
Alpha_rare
dev.off()

