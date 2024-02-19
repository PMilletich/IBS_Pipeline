#####################
#Libraries and Functions 
#####################
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html

#####################
#Load Files 
#####################
setwd("~/Desktop/IBS/")
map = read.csv("IBS_1Year_IBS.csv")


#Recatergorize binned breastfeeding based on Swedish recommendations 
map$Total_Breastfeeding_Binned = ifelse(is.na(map$Total_breastfeeding.mo.), NA, 
                                        ifelse(map$Total_breastfeeding.mo. <= 4, "1-4", 
                                               ifelse(map$Total_breastfeeding.mo. <= 6, "5-6", "7-9")))
map$Exclusive_Breastfeeding_Binned = ifelse(is.na(map$Exclusive_breastfeeding.mo.), NA, 
                                        ifelse(map$Exclusive_breastfeeding.mo. <= 4, "1-4", 
                                               ifelse(map$Exclusive_breastfeeding.mo. <= 6, "5-6", "7-9")))

################################
#Create list of columns to compare
#Remove factors with too many or two few subfactors, and those will category with only two participants 
################################
Confounders_list = colnames(map)

for (current_Confounder in Confounders_list) {
  if (length(unique(map[,current_Confounder])) >= 10 | 
      length(unique(map[,current_Confounder])) == 1) {
    Confounders_list = subset(Confounders_list, Confounders_list != current_Confounder)
  print(current_Confounder)
  } else {
    XXX = data.frame(table(map[,current_Confounder]))

    if (min(XXX$Freq) <= 2) {
      Confounders_list = subset(Confounders_list, Confounders_list != current_Confounder)
      #print(current_Confounder)
    }
  }
}


Confounders_list = subset(Confounders_list, ! Confounders_list %in% 
                            c(     "cdABP", "srABP", "cdsrABP", "cdFG", "IBS",
                                   "Any","abisnr","abp2","abp5","abp8","abp12",
                                   "abp20","k58","k590","k59","k591","k598",
                                   "k599","r10","r104",
                                   "AllAges", "Autoimmune_2_groups", "County", 
                                   "DQ2", "DQ8", "DQ6"))

################################
#Test between controls and all ABP/IBS subgroups 
Case_Confounders = data.frame()

for (current_variable in Confounders_list) {
  set.seed(3)
  #Remove rows with Na's for current variable
  map_subset_current = subset(map, is.na(map[,current_variable]) == F)
  all_current = data.frame()
  if (nrow(data.frame(table(map_subset_current[,current_variable]))) > 1) {
    for (current_group in c("cdABP","srABP", "cdsrABP")) {
      if (current_group == "cdsrABP") {
        subset_df = map_subset_current
        subset_df$IBS = ifelse(subset_df$IBS == "Control", "Control", "cdsrABP")
      } else {
        subset_df = subset(map_subset_current, map_subset_current$IBS %in% c("Control", current_group))
      }
      
      subset_df = subset_df[,c("IBS",current_variable)]
      colnames(subset_df) = c("Group", "Variable")
      current_pvalue = chisq.test(subset_df$Variable, subset_df$Group, 
                              simulate.p.value = T)
      current_pvalue = current_pvalue$p.value
      
      p_row = data.frame("Variable" = current_variable, 
                         "Group" = current_group, 
                         "P" = current_pvalue)
      colnames(p_row) = c("Variable", "Group", "P")
      all_current = rbind(all_current, p_row)
    }
  }
  all_current$padj = p.adjust(all_current$P, "BH" )
  Case_Confounders = rbind(Case_Confounders, all_current)
}

Case_Confounders_Sign = subset(Case_Confounders, Case_Confounders$P <= 0.05)

write.csv(Case_Confounders, "./CSV_Files/ChiSq_all_2.csv")


###################
#Create dataframe of proportions across IBS groups of significant variables
All_data = data.frame()
for (current_column in c("Control", "cdsrABP",
                         "srABP",  "cdABP" )) {
  if (current_column == "cdsrABP") {
    current_data = subset(map, map$IBS != "Control")
  } else {
    current_data = subset(map, map$IBS == current_column)
  }
  print(paste(current_column, nrow(current_data)))
  
  current_col = data.frame()
  for (current_variable in c(unique(Case_Confounders_Sign$Variable)) ){
    current_data_1 = round(table(current_data[,current_variable])/nrow(current_data),3)*100
    current_data_2 = data.frame("Variable" = current_variable,
                                "Current" = current_data_1)
    colnames(current_data_2) = c("Variable", "Subvariable", "Group")
    current_col = rbind(current_col, current_data_2)
  }
  
  while (nrow(current_col) < nrow(All_data)) {
    current_col = rbind(current_col, 
                        data.frame("Variable" = NA, 
                                   "Subvariable" = NA, 
                                   "Group" = NA))
  }
  colnames(current_col) = c("Variable","Subvariable",current_column)
  
  if (nrow(All_data) == 0) {
    All_data = current_col
  } else {
    All_data = cbind(All_data, current_col)
  }
}


