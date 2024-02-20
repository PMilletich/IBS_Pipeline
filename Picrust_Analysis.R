library(ggplot2)
library(phyloseq)

####################################################
#File input
####################################################
setwd("~/Desktop/IBS/")
OG_PS = readRDS("IBS.RDS")
Sample_data = data.frame(sample_data(OG_PS))

#KEGG Data 
Predicted  = read.table("./picrust2-2.4.1/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv", sep = "\t")
colnames(Predicted) = Predicted[1,]
Predicted = Predicted[-1,]
#KEGG Description 
Description = read.csv("./picrust2-2.4.1/picrust2_out_pipeline/KO_metagenome_out/ko_info.csv")

###########################
#Format and normalize Picrust output
###########################
rownames(Predicted) = Predicted$`function`; Predicted$`function` = NULL

Predicted <- mutate_all(Predicted, function(x) as.numeric(as.character(x)))

for (i in 1:ncol(Predicted)) {
  sum_i = sum(Predicted[,i])
  Predicted[,i] = (Predicted[,i]/sum_i)*100
}

Predicted_Samples = data.frame(t(Predicted))
Pathway_Names = colnames(Predicted_Samples)
Predicted_Samples$ID = sub(".", "", rownames(Predicted_Samples))

######################################################
#Find P values across all ABP groups
######################################################
#default values 
current_group ="cdsrABP"
current_pathway = Pathway_Names[1]

p_data = data.frame()
for (current_group in c("cdsrABP","cdABP","srABP")) {
  p_data_1 = data.frame()
  if (current_group =="cdsrABP") {
    IBS_df = subset(Sample_data, Sample_data$IBS != "Control")
    IBS_df$IBS ="cdsrABP"
  } else {
    IBS_df = subset(Sample_data, Sample_data$IBS == current_group)
  }
  
  
  Control = subset(Sample_data, Sample_data$IBS == "Control")
  print(paste(current_group, nrow(Control), nrow(IBS_df)))

  current_Samples = rbind(Control, IBS_df)
  
  Predicted_current = merge(Predicted_Samples, current_Samples, by = "ID")
  for (current_pathway in Pathway_Names) {
    current_data = Predicted_current[,c("ID", "IBS", current_pathway)]
    colnames(current_data)  = c("ID", "IBS", "Pathway")
    current_data$Pathway = as.numeric(current_data$Pathway)
    
    Control = subset(current_data, current_data$IBS == "Control")$Pathway
    IBS = subset(current_data, current_data$IBS != "Control")$Pathway
    
    Control_prev = length(subset(Control, Control>0))/length(Control)
    IBS_prev = length(subset(IBS, IBS>0))/length(IBS)
    
    Control_Mean = mean(Control, na.rm = T)
    IBS_Mean = mean(IBS, na.rm = T)
    
    if (Control_prev >= 0.45 | IBS_prev >= 0.45){
      Wilcox_output = wilcox.test(Pathway ~ IBS, data = current_data)
      W_Pvalue = Wilcox_output$p.value
      
      current_p = data.frame("P" = W_Pvalue, 
                             "Group" = current_group, 
                             "Pathway" = current_pathway,
                             "logDiff" = log(Control_Mean/IBS_Mean))
      p_data_1 = rbind(p_data_1, current_p)
    }
  }
  p_data_1$padj = p.adjust(p_data_1$P, method = "BH")
  p_data = rbind(p_data, p_data_1)
}

######################################################
#Merge for description
set.seed(1)
p_data_2 = merge(p_data, Description, all.x = T)
p_data_sign  = subset(p_data_2, is.na(p_data_2$Description) == F)
p_data_sign  = subset(p_data_sign, p_data_2$Description != "uncharacterized_protein")

#######################################################
#Subset for those sign (p<= 0.05) in all three groups and 0.001 in atleast one
p_data_sign_1 = subset(p_data_sign, p_data_sign$P <= 0.001)
p_data_sign_2 = subset(p_data_sign, p_data_sign$P <= 0.05)

p_data_count = data.frame(table(p_data_sign_2$Pathway))
p_data_count = subset(p_data_count, p_data_count$Freq >= 3)
p_data_count = subset(p_data_count, p_data_count$Var1 %in% p_data_sign_1$Pathway)

p_data_sign_3 = subset(p_data_2, p_data_2$Pathway %in% p_data_count$Var1)
write.csv(p_data_sign_3, "./Images/IBS_4/Picrust_output.csv", row.names = F)

######################################################
#Add spread of relative abundance for all significant predictions
######################################################
graph_data = data.frame()
for (current_group in c("cdsrABP","cdABP","srABP", "Control")) {
  if (current_group =="cdsrABP") {
    IBS_df = subset(Sample_data, Sample_data$IBS != "Control")
    IBS_df$IBS ="cdsrABP"
  } else {
    IBS_df = subset(Sample_data, Sample_data$IBS == current_group)
  }

  current_Samples = IBS_df
  
  Predicted_current = merge(Predicted_Samples, current_Samples, by = "ID")
  for (current_pathway in unique(p_data_sign_3$Pathway)) {
    current_data = Predicted_current[,c("ID", "IBS", current_pathway)]
    colnames(current_data)  = c("ID", "IBS", "Abundance")
    current_data$Abundance = as.numeric(current_data$Abundance)
    current_data$Pathway = current_pathway
    graph_data = rbind(graph_data, current_data)
  }
}

#Reorder X axis 
graph_data$IBS = factor(graph_data$IBS, 
                        levels = c("Control", "cdsrABP", 
                                   "cdABP", "srABP"))
#Create comparison list 
comparison_list = list(c("Control","cdsrABP"), 
                       c("Control","cdABP"), 
                       c("Control","srABP"))

#Merge to have description 
#Posthoc change description names to fit graph better 
graph_data_1 = merge(graph_data, Description)
graph_data_1$Description = gsub("_", " ", graph_data_1$Description)
graph_data_1$Description = gsub(",subfamily,", "\nsubfamily", graph_data_1$Description)
graph_data_1$Description = gsub(",", "\n", graph_data_1$Description)
graph_data_1$Description = gsub(" F ", " F\n", graph_data_1$Description)
graph_data_1$Description = gsub("/formation", "/formation\n", graph_data_1$Description)
graph_data_1$Description = gsub("transport system", 
                                "transport system\n", 
                                graph_data_1$Description)
graph_data_1$Description_1= paste(graph_data_1$Pathway, 
                                  graph_data_1$Description, sep  = "\n")

unique(graph_data_1$Description_1)

#save as image
jpeg("./Images/IBS_4/Picrust_Pathways.jpeg", 
     res = 500, height = 5000, width = 6000)
ggplot(graph_data_1, aes(x = IBS, y = log(Abundance+0.00001), 
                         fill = IBS, color = IBS)) + 
  geom_boxplot()  + theme_bw()+ ylab("log(Relative Abundance)") + 
  xlab(element_blank()) + 
  facet_wrap(~Description_1, scale = "free", ncol = 4) + 
  stat_compare_means(comparisons = comparison_list,
                     vjust = 0.6,step.increase = 0.1,
                     label = "p.signif", hide.ns = T, method = "wilcox")+ 
  theme(legend.position="top") + 
  scale_fill_manual(name = "IBS Group",
                    breaks = c("cdsrABP", "srABP", 
                               "cdABP", "Control"), 
                    values = c("#65ACEA", "#F9D057",
                               "#D05783", "grey25"))+ 
  scale_color_manual(name = "IBS Groups", 
                     breaks = c("cdsrABP", "srABP",
                                "cdABP", "Control"), 
                     values = c("darkblue", "goldenrod",
                                "deeppink4", "grey5"),
                     guide = "none")
dev.off()

