library(phyloseq)
library(microbiomeMarker)
library(metagMisc) #remotes::install_github("vmikk/metagMisc")

#################################
###File Import 
setwd("~/Desktop/IBS/")
physeq_f = readRDS("./IBS.RDS")
physeq_f = tax_glom(physeq_f, "Genus")
sample = read.csv("./IBS_1Year_IBS.csv")

rownames(sample)= sample$ID
sample_data(physeq_f) = sample_data(sample)

Prevalence_threshold = 0.45
#################################################################
#Abundance
#################################################################
ps.RA = transform_sample_counts(physeq_f, function(x) x / sum(x) )

qpcr_data = sample[,c("ID", "copies_16s_per_gram_stool")]
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

#################################################################
#Create Master Dataframe
#################################################################
current_group = "cdsrABP"
lefse_total = data.frame()
aldex_total = data.frame()

for (current_group in c("cdsrABP", "srABP", "cdABP")) {
  current_ps = ps.RA.qpcr
  
  print(current_group)
  #####################################################################################
  #Match Controls to Cases Subset 
  if (current_group != "cdsrABP") {
    IBS_df = subset(sample, sample$IBS == current_group)
  } else {
    IBS_df = subset(sample, sample$IBS != "Control")
    IBS_df$IBS = "cdsrABP"
  }
  Controls_df = subset(sample, sample$IBS == "Control")
  

  IBS_df = rbind(IBS_df, Controls_df)
  rownames(IBS_df) = IBS_df$ID
  print(table(IBS_df$IBS))
  
  #Reorder IBS to ensure same direction for each analysis
  IBS_df$IBS = factor(IBS_df$IBS, 
                      levels = c("Control", current_group))
  sample_data(current_ps) = sample_data(IBS_df)

  current_ps =  phyloseq_rename_with_tax(current_ps, taxrank = "Genus")
  otu_df = data.frame(otu_table(current_ps))
  tax_list = unique(colnames(otu_df))
  otu_df$ID = rownames(otu_df)
  total_df = merge(otu_df, IBS_df)
  
  #Remove genera if prevalence is below threshold 
  #in either abp/ibs group AND controls
  for (current_taxa in tax_list) {
    current_list = total_df[,c(current_taxa, "IBS")]
    Control_prev = subset(current_list, current_list$IBS == "Control")
    Control_prev = length(subset(Control_prev[,current_taxa], 
                                 Control_prev[,current_taxa] > 0))/length(Control_prev[,current_taxa])
    ABP_prev = subset(current_list, current_list$IBS != "Control")
    ABP_prev = length(subset(ABP_prev[,current_taxa], 
                                 ABP_prev[,current_taxa] > 0))/length(ABP_prev[,current_taxa])
    
    if (Control_prev < Prevalence_threshold & ABP_prev < Prevalence_threshold) {
      otu_df[,current_taxa] = NULL
    }
  }
  otu_df$ID = NULL
  otu_table(current_ps) = otu_table(otu_df, taxa_are_rows = F)

  #########################################################
  #LEFSE
  #########################################################
  mm_lefse <- run_lefse(
    current_ps,
    taxa_rank = "Genus",
    group = "IBS",
    norm = "none"
  )
  
  lefse_df = data.frame(marker_table(mm_lefse))
  lefse_df$Group = current_group
  lefse_total = rbind(lefse_total, lefse_df)
  
  #########################################################
  #ALDEX
  #########################################################
  mm_aldex = suppressMessages(run_aldex(
    current_ps,
    "IBS", taxa_rank = "Genus",
    transform = c("identity"), norm = "none",
    method = c("t.test"),
    pvalue_cutoff = 0.05,
    mc_samples = 128,
    denom = c("all"),
    paired = FALSE
  ))
  
  aldex_df = data.frame(marker_table(mm_aldex))
  aldex_df$Group = current_group
  aldex_total = rbind(aldex_total, aldex_df)
  
}


write.csv(lefse_total, "./Images/IBS_4/Lefse.csv", row.names = F)
write.csv(aldex_total, "./Images/IBS_4/Aldex2.csv", row.names = F)

table(lefse_total$feature)
lefse_total$ef_lda_1 = ifelse(lefse_total$enrich_group == "Control", 
                            lefse_total$ef_lda*-1, lefse_total$ef_lda)

lefse_plot = ggplot(lefse_total, aes(x = ef_lda_1, y = feature, 
                                     color = Group, 
                                     size = -log(padj))) + 
  geom_point(alpha = 0.75) + 
    theme_bw()+ 
  xlab("ef_lda") + 
  geom_vline(xintercept = 0) + 
  scale_fill_manual(breaks = c("cdsrABP", "srABP", 
                               "cdABP", "Control"), 
                    values = c("#65ACEA", "#F9D057",
                               "#D05783", "grey25"))  + 
    scale_color_manual(breaks = c("cdsrABP", "srABP", 
                                  "cdABP", "Control"), 
                       values = c("#65ACEA", "#F9D057",
                                  "#D05783", "grey25")); lefse_plot

table(aldex_total$feature)
# aldex_total$ef_aldex = ifelse(aldex_total$enrich_group == "Control", 
#                             aldex_total$ef_aldex, aldex_total$ef_aldex*-1)

Aldex_plot = ggplot(aldex_total, aes(x = ef_aldex, y = feature, color = Group, 
                        size = -log(padj))) + 
  geom_point(alpha = 0.75) + 
  theme_bw()+ 
  geom_vline(xintercept = 0) + 
  scale_fill_manual(breaks = c("cdsrABP", "srABP", 
                               "cdABP", "Control"), 
                    values = c("#65ACEA", "#F9D057",
                               "#D05783", "grey25"))  + 
  scale_color_manual(breaks = c("cdsrABP", "srABP", 
                                "cdABP", "Control"), 
                     values = c("#65ACEA", "#F9D057",
                                "#D05783", "grey25")) 

jpeg("./Images/IBS_4/Differential_Abundance.jpeg", res = 500, 
     height = 2500, width = 4000)
ggarrange(lefse_plot, Aldex_plot,
          ncol = 2, common.legend = T)
dev.off()

