#####################
#Libraries and Functions 
#####################
library(vegan) #https://cran.r-project.org/web/packages/vegan/index.html
library(plyr) #http://myweb.facstaff.wwu.edu/minerb2/biometrics/plyr.html
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
'%notin%' = Negate('%in%')#https://stackoverflow.com/questions/38351820/negation-of-in-in-r


#https://www.webpages.uidaho.edu/~stevel/519/How%20to%20get%20correlation%20between%20two%20categorical%20variable%20and%20a%20categorical%20variable%20and%20continuous%20variable.html
#https://rpubs.com/hoanganhngo610/558925
#####################
#Load Files 
#####################
#Sample Data and Phyloseq from Previous Subsetting Step 
map = read.csv("IBS_1Year.csv")

################################
#Create list of columns to compare 
################################
Confounders_list = colnames(map)

for (current_Confounder in Confounders_list) {
  if (length(unique(map[,current_Confounder])) >= 10 | length(unique(map[,current_Confounder])) == 1) {
    Confounders_list = subset(Confounders_list, Confounders_list != current_Confounder)
    #print(current_Confounder)
  } else {
    XXX = data.frame(table(map[,current_Confounder]))
    
    if (min(XXX$Freq) <= 4) {
      Confounders_list = subset(Confounders_list, Confounders_list != current_Confounder)
      #print(current_Confounder)
    }
  }
}

Confounders_list = subset(Confounders_list, ! Confounders_list %in% 
                            c("abp_Multiple","k_Combined","r_Combined",
                              "Any","abisnr","abp2","abp5","abp8","abp12",
                              "abp20","k58","k590","k59","k591","k598",
                              "k599","r10","r104",
                              "County", "AllAges",
                              "Exclusive_Breastfeeding_Binned_Old" ,"Total_Breastfeeding_Binned_Old",
                              "Intro_formula_Binned_Old","Intro_CowsMilk_Binned_Old",
                              "Intro_Gluten_Binned_Old","Intro_Fish_Binned_Old" ))

################################
#### Factors impacting Autoimmune Development 
################################
#create new dataframes for Cases or FAID_Nos

environmental_list = Confounders_list

current_variable = "Siblings_at_birth"
Case_Confounders = data.frame()

map$abp_Combined = ifelse(map$abp_Multiple == "No", "No", "Any")
table(map$Any)

for (current_variable in environmental_list) {
  set.seed(3)
  #Create a new dataframe of subjects with information on the current variable 
  map_subset_current = subset(map, is.na(map[,current_variable]) == F)
  map_subset_current$abp_Combined = ifelse(map_subset_current$abp_Multiple == "No", "No", "Any")
  if (nrow(data.frame(table(map_subset_current[,current_variable]))) > 1) {
    None = subset(map_subset_current, map_subset_current$Any == "No")[,current_variable]
    Any = subset(map_subset_current, map_subset_current$Any == "Yes")[,current_variable]
    R = subset(map_subset_current, map_subset_current$r_Combined == "Any")[,current_variable]
    K = subset(map_subset_current, map_subset_current$k_Combined == "Any")[,current_variable]
    abp = subset(map_subset_current, map_subset_current$abp_Combined == "Any")[,current_variable]
    
    #Any
    Any_pvalue = data.frame("Group" = c(rep("None", length(None)),
                                       rep("Any", length(Any))), 
                           "Variable" = c(None, Any))
    Any_pvalue = chisq.test(Any_pvalue$Variable, Any_pvalue$Group, 
                            simulate.p.value = T)
    Any_pvalue = Any_pvalue$p.value
    
    #r_Combined
    set.seed(1)
    R_pvalue = data.frame("Group" = c(rep("None", length(None)),
                                        rep("R", length(R))), 
                            "Variable" = c(None, R))
    R_pvalue = chisq.test(R_pvalue$Variable, R_pvalue$Group, 
                            simulate.p.value = T)
    R_pvalue = R_pvalue$p.value
    
    #abp
    abp_pvalue = data.frame("Group" = c(rep("None", length(None)),
                                        rep("abp", length(abp))), 
                            "Variable" = c(None, abp))
    abp_pvalue = chisq.test(abp_pvalue$Variable, abp_pvalue$Group, 
                            simulate.p.value = T)
    abp_pvalue = abp_pvalue$p.value
    
    #Any
    K_pvalue = data.frame("Group" = c(rep("None", length(None)),
                                        rep("K", length(K))), 
                            "Variable" = c(None, K))
    K_pvalue = chisq.test(K_pvalue$Variable, K_pvalue$Group, 
                            simulate.p.value = T)
    K_pvalue = K_pvalue$p.value
    
    #Bind the pvalue and variable to master dataframe 
    current_row = data.frame("Variable" = current_variable, 
                             "Any" = Any_pvalue,
                             "r_Combined" = R_pvalue,
                             "k_Combined" = K_pvalue, 
                             "abp_Multiple" = abp_pvalue)
    Case_Confounders = rbind(Case_Confounders, current_row)
  }
}

Case_Confounders$Any_FDR = p.adjust(Case_Confounders$Any, method = "BH")
Case_Confounders$abp_FDR = p.adjust(Case_Confounders$abp_Multiple, method = "BH")
Case_Confounders$k_FDR = p.adjust(Case_Confounders$k_Combined, method = "BH")
Case_Confounders$r_FDR = p.adjust(Case_Confounders$r_Combined, method = "BH")

Case_Confounders_Sign = subset(Case_Confounders, Case_Confounders$Any <= 0.05 | 
                                 Case_Confounders$abp_Multiple <= 0.05 | 
                                 Case_Confounders$k_Combined <= 0.05 | 
                                 Case_Confounders$r_Combined <= 0.05)

#write.csv(Case_Confounders_Sign, "./CSV_Files/ChiSq_1.csv")
