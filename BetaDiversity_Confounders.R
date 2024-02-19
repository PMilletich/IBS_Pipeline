#####################
#Libraries and Functions 
#####################
library(vegan) #https://cran.r-project.org/web/packages/vegan/index.html
library(plyr) #http://myweb.facstaff.wwu.edu/minerb2/biometrics/plyr.html
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

#####################
#Load Files 
#####################
setwd("~/Desktop/IBS/")
OG_PS = readRDS("IBS.RDS")
OG_PS = tax_glom(OG_PS, 'Genus')
map = read.csv("IBS_1Year_IBS.csv")

################################
#Create list of columns to compare 
#Remove columns with too many variables or low coints
################################
Confounders_list = colnames(map)

for (current_Confounder in Confounders_list) {
  if (length(unique(map[,current_Confounder])) >= 10 | length(unique(map[,current_Confounder])) == 1) {
    Confounders_list = subset(Confounders_list, Confounders_list != current_Confounder)
  } else {
    XXX = data.frame(table(map[,current_Confounder]))
    
    if (min(XXX$Freq) <= 4) {
      Confounders_list = subset(Confounders_list, Confounders_list != current_Confounder)
    }
  }
}

Confounders_list = subset(Confounders_list, ! Confounders_list %in% 
                            c(     "cdABP", "srABP", "cdsrABP", "cdFG", "IBS",
                                   "Any","abisnr","ID",
                                   "abp2","abp5","abp8","abp12","abp20",
                                   "r10","r104",
                                   "AllAges", "Autoimmune_2_groups", "County", 
                                   "DQ2", "DQ8", "DQ6"))

#################################################################
#Relative Abundance 
#################################################################
ps.RA = transform_sample_counts(OG_PS, function(x) x / sum(x) )

#################################################################
#Absolute Abundance 
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

################################
#### Factors impacting Binomial Beta Diversity
#### Takes a long time to run, depending on number of columns comparing 
################################
current_ps = ps.qpcr 
rownames(map) = map$ID
sample_data(current_ps) = sample_data(map)
Case_Confounders = data.frame()
current_variable = "Age.collected.month"
for (current_variable in Confounders_list) {
  set.seed(5)
  ps_Variable = current_ps
  #Create a new sample dataframe and remove subjects lacking the current_variable
  metadata_Variable = as(sample_data(ps_Variable), "data.frame")
  metadata_Variable = subset(metadata_Variable, is.na(metadata_Variable[,current_variable]) == F)
  metadata_Variable[,current_variable] = lapply(metadata_Variable[current_variable], as.character)
  table(metadata_Variable[,current_variable])
  if (nrow(data.frame(table(metadata_Variable[,current_variable]))) > 1) {
    #Load in subsetted data to PS
    sample_data(ps_Variable) = sample_data(metadata_Variable)
    
    #Create a distance matrix from the phyloseq object 
    #https://joey711.github.io/phyloseq/distance.html
    dist.uf <- phyloseq::distance(ps_Variable, method = "binomial",na.rm = TRUE)
    
    #Run anova between the distance matrix and the current_variable
    #https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis
    ANOVA = adonis2(as.formula(paste("dist.uf ~ ", current_variable, sep = "")),
                    data = metadata_Variable)
    
    #Extract the Pvalue
    Pvalue = ANOVA$`Pr(>F)`[1]
    R2 = round(ANOVA$R2[1],3)
    DF = ANOVA$Df[1]
    
    #Append the Pvalue and the variable to the master dataframe 
    current_row = data.frame("Variable" = current_variable,
                             "DF" = DF, 
                             "R2" = R2, 
                             "Pvalue" = Pvalue)
    print(current_row)
    Case_Confounders = rbind(Case_Confounders, current_row)
  }
}

Case_Confounders$padj = p.adjust(Case_Confounders$Pvalue, "BH")

#Save dataframe to csv 
write.csv(Case_Confounders, "./CSV_Files/Beta_Diversity_3.csv", row.names = F)
