library(ggplot2)
library(phyloseq)
library(dplyr)
library(ggpubr)

####################################################
#File input
####################################################
setwd("~/Desktop/IBS/")
OG_PS = readRDS("IBS.RDS")
Sample_data = data.frame(sample_data(OG_PS))

#KEGG Data 
Predicted  = read.table("./picrust2-2.4.1/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv", sep = "\t")
colnames(Predicted) = Predicted[1,]
Predicted = Predicted[-1,]
#Metacyc Description 
#https://github.com/picrust/picrust2/tree/master/picrust2/default_files/description_mapfiles
Description = read.csv("./picrust2-2.4.1/picrust2_out_pipeline/pathways_out/metacyc_pathways_info.csv")

###########################
#Format and normalize Picrust output
###########################
rownames(Predicted) = Predicted$pathway; Predicted$pathway = NULL

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
p_data_sign  = subset(p_data_sign, grepl("uncharacterized_protein", p_data_2$Description) == F)

#######################################################
#Subset for those sign (p<= 0.05) in all three groups
p_data_sign_2.1 = subset(p_data_sign, p_data_sign$P < 0.05)
length(unique(p_data_sign_2.1$Pathway))
p_data_sign_2 = subset(p_data_sign, p_data_sign$P < 0.1)

p_data_count = data.frame(table(p_data_sign_2$Pathway))
p_data_count = subset(p_data_count, p_data_count$Freq >= 3)

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

graph_data_1$Description_1= paste(graph_data_1$Pathway,
                                  graph_data_1$Description, 
                                  sep  = "\n")

names = unique(graph_data_1$Description_1)
names
#save as image
jpeg("./Images/IBS_4/Picrust_Pathways.jpeg", 
     res = 500, height = 5000, width = 6000)
ggplot(graph_data_1, aes(x = IBS, y = log(Abundance+0.00001), 
                         fill = IBS, color = IBS)) + 
  geom_boxplot()  + theme_bw()+ ylab("log(Relative Abundance)") + 
  xlab(element_blank()) + 
  facet_wrap(~Description_1, scale = "free", ncol = 3) + 
  stat_compare_means(comparisons = comparison_list, 
                     hide.ns = T, method = "wilcox")+ 
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

############################################################################
#Confounders 
Total_data = merge(Sample_data, Predicted_Samples, by = "ID")
Heterogeneous = c("Sex", "DR15.DQ602", "DR4.DQ8")
Compare_list= unique(graph_data_1$Pathway)
Sign_data = data.frame()
for (current_variable in Heterogeneous) {
  pvalue_data = data.frame()
  current_data = Total_data 
  current_data = subset(current_data, is.na(current_data[,current_variable]) == F)
  current_data[,current_variable] = as.factor(current_data[,current_variable])
  
  current_pathway= Compare_list[1]
  #If there is more than one sub-category 
  if (length(unique(current_data[,current_variable]))>1) {
    for (current_pathway in Compare_list) {
      current_pathway = gsub("/", ".", current_pathway)
      #Determine normalcy
      norm = shapiro.test(current_data[,current_pathway])
      if (norm$p.value > 0.05) { #normal
        ANOVA = aov(as.formula(paste(current_pathway, "~", current_variable)),
                    data = current_data)
        Pvalue = summary(ANOVA)[[1]][["Pr(>F)"]][1]
        
      } else {
        ANOVA = kruskal.test(as.formula(paste(current_pathway, "~", current_variable)),
                             data = current_data)
        Pvalue = ANOVA$p.value
      }
      
      #Run kruskal.test, extract and save p-value
      P_row = data.frame("Pathway" = current_pathway, 
                         "Variable" = current_variable,  
                         "Pvalue" =  Pvalue)
      pvalue_data = rbind(pvalue_data, P_row)
    }
  }
  Sign_data = rbind(Sign_data, pvalue_data)
}
  
Sign_data$FDR = p.adjust(Sign_data$Pvalue, "fdr")
Sign_data$FDR_1 = p.adjust(Sign_data$Pvalue, method = "fdr")
Sign_data = subset(Sign_data, Sign_data$Pvalue <= 0.05)

#####################################
#Graphing  
table(Sign_data$Variable); table(Sign_data$Pathway)
current_variable = "Sex"
current_Pathway = "PWY.6629"
for (current_variable in unique(Sign_data$Variable)) {
  graphing_data = data.frame()
  
  current_subset = subset(Sign_data,Sign_data$Variable == current_variable)
  for (current_Pathway in unique(current_subset$Pathway)) {
    current_count = Total_data
    current_subset_1 = current_count[,c(current_variable, current_Pathway)]
    current_subset_1$Factor = current_variable
    current_subset_1$Pathway = current_Pathway
    current_subset_1$Group = "All"
    colnames(current_subset_1) = c("Subfactor", "Value", "Factor", "Pathway", "Group")
    current_subset_1= subset(current_subset_1, is.na(current_subset_1$Subfactor) == F)
    graphing_data = rbind(graphing_data, current_subset_1)
  }
  graphing_data$Pathway = gsub("_", " ", graphing_data$Pathway)
  graphing_data$Pathway = gsub("Escherichia.Shigella", "Escherichia/Shigella", graphing_data$Pathway)
  
  
  col_list = c("red4","orangered","orange","gold","lightgreen","chartreuse2",
               "darkgreen","turquoise3", "cyan1","dodgerblue","blue", "darkblue",
               "purple2", "magenta4","hotpink2","orchid1", "plum1", 
               "rosybrown1", "tan","chocolate", "tan4", "slategray2", "azure4")
  #col_list = rainbow(length(Compare_list))
  set.seed(11)
  col_list = sample(col_list)

  c_plot = suppressWarnings(ggplot(graphing_data, aes(x = Subfactor,
                                                      Group = Subfactor, 
                                                      y = log2(Value+1))) + 
                              ylab("log2(Relative Abundance)") + 
                              geom_boxplot(color = "black", data = graphing_data, 
                                           aes(fill = Pathway)) + theme_bw() + 
                              # geom_point(position = position_jitterdodge(jitter.width = 0.8),
                              #            alpha = 0.15, size = 1) + 
                              scale_fill_manual(breaks = Compare_list, 
                                                values = col_list) + 
                              ggtitle(current_variable_1) + 
                              facet_wrap(~Pathway, scale = "free", nrow = 1) + 
                              theme(axis.title.x = element_blank(),
                                    legend.position = "none") + 
                              stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), 
                                                 size = 5, vjust = 1))
  suppressWarnings(print(c_plot))
  assign(paste(current_variable, "_plot", sep = ""), c_plot)
  print(paste(current_variable, "_plot", sep = ""))
}

jpeg("./Images/IBS_4/Environmental_Pathway.jpeg", res = 500, 
     height = 2000, width = 4000)
Sex_plot
dev.off()
