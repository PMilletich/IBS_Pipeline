#Create a long dataframe of the IBS data 
#####################
#Load Files 
#####################
library(ggplot2) #https://ggplot2.tidyverse.org/
library(ggpubr) #https://cran.r-project.org/web/packages/ggpubr/index.html
library(phyloseq) #https://joey711.github.io/phyloseq/

setwd("~/Desktop/IBS/")
OG_PS = readRDS("IBS.RDS")
OG_PS = tax_glom(OG_PS, "Genus")
OG_taxa = data.frame(tax_table(OG_PS))
map = read.csv("IBS_1Year_IBS.csv")
table(map$IBS)

set.seed(1)
#################################################################
#Relative Abundance 
ps.RA = transform_sample_counts(OG_PS, function(x) x / sum(x) )

#################################################################
#Absolute Abundance 
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

####################################################
#Create a sample data and phyloseq object all containing "cdsrABP"
####################################################
#OG data with Controls, cdABP, srABP
long_sample = map
rownames(long_sample) = long_sample$ID
long_otu = data.frame(otu_table(ps.qpcr))

##################################################
#Add cdsrABP
current_subset = subset(map, map$IBS != "Control" )
current_subset$IBS = "cdsrABP"
rownames(current_subset) = current_subset$ID
#Create phyloseq object with new sample data; [ 224 taxa and 208 samples ]
current_PS = ps.qpcr
sample_data(current_PS) = sample_data(current_subset)

#Append the Sample data
rownames(current_subset) = paste(current_subset$ID, "cdsrABP", sep = "_")
long_sample = rbind(long_sample, current_subset)

#Create OTU and append to total otu
current_otu_data = data.frame(otu_table(current_PS))
rownames(current_otu_data) = paste(rownames(current_otu_data), "cdsrABP", sep = "_")
long_otu = rbind(long_otu, current_otu_data)
  
#Create new phyloseq object; [ 224 taxa and 684 samples ]
long_sample$IBS = factor(long_sample$IBS, 
                         levels = c("Control", "cdsrABP", 'cdABP', 'srABP'))
long_ps = phyloseq(otu_table(long_otu, taxa_are_rows = F), 
                   tax_table(as.matrix(OG_taxa)), 
                   sample_data(long_sample))

########################################################################
#Rarefy 
########################################################################
#Determine threshold for rarefying 
samplesums = data.frame(sample_sums(long_ps))
samplesums$ID =rownames(samplesums)
summary(samplesums[,1])
samplesums_1 = subset(samplesums, samplesums[,1] <= 200000)

Physeq_rare = rarefy_even_depth(long_ps, sample.size = 200000,
                                rngseed = 13, replace = F, 
                                trimOTUs = F, verbose = T)

################################################################
#Original 
################################################################
long_sample_OG_ps = Physeq_rare

my_comparisons_OG <- list( c("Control", "cdsrABP"),
                           c("Control", "cdABP"),
                           c("Control", "srABP"),
                           c("srABP", "cdABP"),
                           c("cdsrABP", "cdABP"),
                           c("cdsrABP", "srABP"))


Alpha_rare_OG = plot_richness(long_sample_OG_ps, x="IBS", measures=c("Observed", "Shannon"),
                           color = "IBS") + 
  stat_compare_means(comparisons = my_comparisons_OG, 
                     #method = "wilcox", 
                     step.increase = 0.1, tip.length = 0.01) +  
  scale_color_manual(name = "IBS Groups",
                     breaks = c("cdsrABP", "srABP", "cdABP", "Control"), 
                     values = c("#65ACEA", "#F9D057", "#D05783", "grey25")) +
  geom_jitter(alpha = 0.5, width = 0.4, height = 0) +
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        axis.title.x = element_blank(), 
        legend.position = "top",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
Alpha_rare_OG$layers = Alpha_rare_OG$layers[-1];Alpha_rare_OG

jpeg("./Images/IBS_4/Alpha_rare_Subset_Genus_1.jpeg", res = 600,
     width = 4000, height = 3000)
Alpha_rare_OG
dev.off()

