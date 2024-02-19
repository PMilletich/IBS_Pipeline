#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(pime)
library(ggplot2)
library(ggpubr) #ggarrange
library(reshape2)

#################################
###File Import 
setwd("~/Desktop/IBS/")

set.seed(1)

physeq_f = readRDS("./IBS.RDS")
physeq_f = tax_glom(physeq_f, "Genus")
sample = read.csv("./IBS_1Year_IBS.csv")

rownames(sample)= sample$ID

sample_data(physeq_f) = sample_data(sample)
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


##############################################
#Run through all IBS/ABP Groups
##############################################
table(sample$IBS)
current_group= "srABP"
All_MDA = data.frame()

for (current_group in c("srABP", "cdABP", "cdsrABP")) {
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

  PIME_PS = ps.RA.qpcr
  sample_data(PIME_PS) = sample_data(IBS_df)
  set.seed(1)
  
  ##############################################
  #Original OOB with no filtering
  print(pime.oob.error(PIME_PS, "IBS"))
  ##############################################
  #Run iterations to determine best filtering prevalence
  per_variable_obj_IBS= pime.split.by.variable(PIME_PS, "IBS")
  prevalences_IBS=pime.prevalence(per_variable_obj_IBS)
  rep.oob.error= pime.oob.replicate(prevalences_IBS, 
                                    "IBS",  bootstrap = 50, parallel = TRUE)
  
  #create dataframe of oob at all thresholds for all iterations
  rep.oob.df = data.frame(t(rep.oob.error$`Results table`))
  rep.oob.df$Variable = rownames(rep.oob.df)
  #Turn to long dataframe using prevalence and iterations
  rep.oob.df = melt(rep.oob.df, id.vars = "Variable")
  
  #posthoc add color for selected threshold of interest 
  rep.oob.df$Importance = ifelse(rep.oob.df$Variable == "Prevalence45%",
                                 current_group, "NA")
  #remove extra characters from xaxis variables 
  rep.oob.df$Variable = gsub("Prevalence", "", 
                             rep.oob.df$Variable)
  rep.oob.df$Variable = gsub("%", "", 
                             rep.oob.df$Variable)
  #reorder x axis variables to make sure it's in numerical order 
  rep.oob.df$Variable = factor(rep.oob.df$Variable,
                               levels = c("5","10", "15", "20", "25",
                                          "30", "35", "40", "45", "50",
                                          "55", "60", "65", "70", "75",
                                          "80", "85", "90", "95"))
  rep.oob.df$value= rep.oob.df$value*100
  #create boxplots of OOB for all iterations and prevalences 
  oob_plot = ggplot(rep.oob.df, aes(x = Variable, y = value,
                                    color = Importance)) + 
    geom_hline(yintercept = 5) + 
    geom_boxplot() + theme_bw() + ylim(0,50) + 
    scale_color_manual(breaks = c("cdsrABP", "srABP", 
                                  "cdABP", "Control"), 
                       values = c("#65ACEA", "#F9D057",
                                  "#D05783", "grey25"),
                       drop = F) + 
    ggtitle(current_group) + 
    xlab("Prevalence Threshold") + ylab("OOB")+ 
    theme(legend.position = "none")
  
  plot(oob_plot)
  #Assign plot to IBS group variable 
  assign(paste(current_group, "_oob", sep = ""), oob_plot)
  
  
  #Use random iteration for MDA
  set.seed(42)
  best.prev_IBS=pime.best.prevalence(prevalences_IBS, "IBS")
  #use posthoc threshold that matches line 87
  current_imp_IBS=best.prev_IBS$`Importance`$`Prevalence 45`

  
  #Create dataframe of MDA genera 
    current_imp_IBS$Group = current_group
    current_imp_IBS = current_imp_IBS[,c(current_group, "Control",
                                         "MeanDecreaseAccuracy", 
                                         "Genus", "Group")]
    colnames(current_imp_IBS) = c("Yes", "No", "MeanDecreaseAccuracy", "Genus", "Group")
    All_MDA = rbind(All_MDA, current_imp_IBS)
}

#Save All boxplots 
jpeg("./Images/IBS_4/PCOA_4Groups_4_OOB.jpeg", 
     res = 500, height = 3000, width = 4000)
ggarrange(cdsrABP_oob, 
          ggarrange(srABP_oob, cdABP_oob, ncol = 2), 
          nrow = 2)
dev.off()

#subset based on mda score
All_MDA_Sign_1 = subset(All_MDA, All_MDA$MeanDecreaseAccuracy > 0)
write.csv(All_MDA_Sign_1, "./Images/IBS_4/PIME_All_4_1.csv")

#how many times the genera are in all three groups
df = data.frame(table(All_MDA_Sign_1$Genus)) 


