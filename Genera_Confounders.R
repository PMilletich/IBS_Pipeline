#####################
#Libraries and Functions 
#####################
library(vegan) #https://cran.r-project.org/web/packages/vegan/index.html
library(plyr) #http://myweb.facstaff.wwu.edu/minerb2/biometrics/plyr.html
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
library(ggplot2)
library(ggpubr)
library(metagMisc)

#####################
#Load Files 
#####################
setwd("~/Desktop/IBS/")
PS = readRDS("IBS.RDS")
PS = tax_glom(PS, "Genus")
OG_taxa = data.frame(tax_table(PS))

sample_data = data.frame(sample_data(PS))

###################
#Significant Taxa
###################
PIME = read.csv("./Images/IBS_4/PIME_All_4_1.csv")
Aldex = read.csv("./Images/IBS_4/Aldex2.csv")
Lefse = read.csv("./Images/IBS_4/Lefse.csv")

PIME_Lefse = intersect(PIME$Genus, Lefse$feature)
PIME_Aldex = intersect(PIME$Genus, Aldex$feature)
Lefse_Aldex = intersect(Lefse$feature, Aldex$feature)
Compare_list = unique(c(PIME_Aldex, PIME_Lefse))
Compare_list = unique(c(Compare_list, Lefse_Aldex))

################################
#Create list of columns to compare 
################################
Heterogeneous = read.csv("./CSV_Files/ChiSq_all_2.csv")
Heterogeneous = subset(Heterogeneous, Heterogeneous$P <= 0.05)
Heterogeneous = c(unique(Heterogeneous$Variable))

#################################################################
#Relative Abundance 
#################################################################
ps.RA = transform_sample_counts(PS, function(x) x / sum(x) )

#################################################################
#qPCR 
#################################################################
qpcr_data = sample_data[,c("ID", "copies_16s_per_gram_stool")]
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

#################################################################
#Create Master Dataframe
#################################################################
final_data = data.frame()
current_ps = ps.qpcr

current_Taxa = data.frame(tax_table(current_ps))
current_Taxa = subset(current_Taxa, current_Taxa$Genus %in% Compare_list)
tax_table(current_ps) = tax_table(as.matrix(current_Taxa))

current_ps = phyloseq_rename_with_tax(current_ps, "Genus")
count_table = data.frame(otu_table(current_ps))
count_table$ID = rownames(count_table)

#Merge OTU with Sample
Total_data = merge(count_table, sample_data)

############################################
#Statistical Testing
############################################
current_variable = "Sex"
#current_variable = Confounders_list[1]
Sign_data = data.frame()
#Test significance of all metadata on Genera 
for (current_variable in Heterogeneous) {
  pvalue_data = data.frame()
  current_data = Total_data 
  #Remove NAs and turn into factors 
  current_data = subset(current_data, is.na(current_data[,current_variable]) == F)
  current_data[,current_variable] = as.factor(current_data[,current_variable])
  
  for (current_genus in Compare_list) {
    current_genus = gsub("/", ".", current_genus)
    #Determine normalcy
    norm = shapiro.test(current_data[,current_genus])
    if (norm$p.value > 0.05) { #normal
      ANOVA = aov(as.formula(paste(current_genus, "~", current_variable)),
                  data = current_data)
      Pvalue = summary(ANOVA)[[1]][["Pr(>F)"]][1]
      
    } else {
      ANOVA = kruskal.test(as.formula(paste(current_genus, "~", current_variable)),
                           data = current_data)
      Pvalue = ANOVA$p.value
    }
    
    #Run kruskal.test, extract and save p-value
    P_row = data.frame("Genus" = current_genus, 
                       "Variable" = current_variable,  
                       "Pvalue" =  Pvalue)
    pvalue_data = rbind(pvalue_data, P_row)
  }
  Sign_data = rbind(Sign_data, pvalue_data)
}

Sign_data$FDR_1 = p.adjust(Sign_data$Pvalue, method = "fdr")
Sign_data = subset(Sign_data, Sign_data$Pvalue <= 0.05)

#####################################
#Create long data of absolute abundance for significant factors   
for (current_variable in unique(Sign_data$Variable)) {
  graphing_data = data.frame()
  
  current_subset = subset(Sign_data,Sign_data$Variable == current_variable)
  for (current_genus in unique(current_subset$Genus)) {
    current_count = Total_data
    current_subset_1 = current_count[,c(current_variable, current_genus)]
    current_subset_1$Factor = current_variable
    current_subset_1$Genus = current_genus
    current_subset_1$Group = "All"
    colnames(current_subset_1) = c("Subfactor", "Value", "Factor", "Genus", "Group")
    current_subset_1= subset(current_subset_1, is.na(current_subset_1$Subfactor) == F)
    graphing_data = rbind(graphing_data, current_subset_1)
  }
  graphing_data$Genus = gsub("_", " ", graphing_data$Genus)
  graphing_data$Genus = gsub("Escherichia.Shigella", "Escherichia/Shigella", graphing_data$Genus)
  
  
  #Create unique colors that would match genera across variour factors 
  col_list = c("red4","orangered","orange","gold","lightgreen","chartreuse2",
               "darkgreen","turquoise3", "cyan1","dodgerblue","blue", "darkblue",
               "purple2", "magenta4","hotpink2","orchid1", "plum1", 
               "rosybrown1", "tan","chocolate", "tan4", "slategray2", "azure4")
  set.seed(3)
  col_list = sample(col_list)
  
  c_plot = ggplot(graphing_data, aes(x = Subfactor, Group = Subfactor, 
                                     y = log2(Value+1))) + 
    ylab("log2(Abundance)") + 
    geom_boxplot(color = "black", data = graphing_data, 
                 aes(fill = Genus)) + theme_bw() + 
    scale_fill_manual(breaks = Compare_list, 
                      values = col_list) + 
    ggtitle(current_variable_1) + 
    facet_wrap(~Genus, scale = "free", nrow = 1) + 
    theme(axis.title.x = element_blank(),
          legend.position = "none") + 
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), 
                       size = 5, vjust = 1)
  assign(paste(current_variable, "_plot", sep = ""), c_plot)
  print(paste(current_variable, "_plot", sep = ""))
}

jpeg("./Images/IBS_4/Environmental_Genus.jpeg", res = 500, 
     height = 2000, width = 4000)
Sex_plot
dev.off()

