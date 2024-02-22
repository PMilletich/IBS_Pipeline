#####################
#Libraries and Functions 
#####################
#library(vegan) #https://cran.r-project.org/web/packages/vegan/index.html
#library(plyr) #http://myweb.facstaff.wwu.edu/minerb2/biometrics/plyr.html
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
library(epitools)
library(ggplot2)
library(scales)

#####################
#Load Files 
#####################
setwd("~/Desktop/IBS/")
map = read.csv("IBS_1Year_IBS.csv")


#Recatergorize binned breastfeeding based on swedish recommendations 
map$Total_Breastfeeding_Binned = ifelse(is.na(map$Total_breastfeeding.mo.), NA, 
                                        ifelse(map$Total_breastfeeding.mo. <= 4, "1-4", 
                                               ifelse(map$Total_breastfeeding.mo. <= 6, "5-6", "7-9")))
map$Exclusive_Breastfeeding_Binned = ifelse(is.na(map$Exclusive_breastfeeding.mo.), NA, 
                                        ifelse(map$Exclusive_breastfeeding.mo. <= 4, "1-4", 
                                               ifelse(map$Exclusive_breastfeeding.mo. <= 6, "5-6", "7-9")))

################################
#Create list of columns to compare 
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

#### Factors impacting Autoimmune Development 
################################
#create new dataframes for Cases or FAID_Nos

Confounders_list = Confounders_list

current_variable = "Siblings_at_birth"
Case_Confounders = data.frame()

for (current_variable in Confounders_list) {
  set.seed(3)
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

Case_Confounders_2 = subset(Case_Confounders, Case_Confounders$Variable %in% 
                              Case_Confounders_Sign$Variable)
write.csv(Case_Confounders, "./CSV_Files/ChiSq_all_2.csv")

All_data = data.frame()
for (current_column in c("cdsrABP",
                         "srABP",  "cdABP" )) {
  if (current_column == "cdsrABP") {
    current_data = map
    current_data$IBS = ifelse(current_data$IBS == "Control", 
                              "Control", "cdsrABP")
  } else {
    current_data = subset(map, map$IBS %in% c(current_column, "Control"))
  }

  current_data$IBS = factor(current_data$IBS, 
                            levels= c(current_column, "Control"))
  
  current_col = data.frame()
  for (current_variable in c(unique(Case_Confounders_Sign$Variable)) ){
   print(current_variable)
     current_table = table(current_data$IBS, current_data[,current_variable])
    OR = oddsratio(current_table)
    print(OR$data)
    p = OR$p.value[2,3]
    OR_value = OR$measure[2,1]
    UL_value = OR$measure[2,3]
    LL_value = OR$measure[2,2]
    
    current_data_2 = data.frame("Group" = current_column,
                                "Variable" = current_variable,
                                "OR" = OR_value, 
                                "UL" = UL_value, 
                                "LL" = LL_value
                                )

    All_data = rbind(All_data, current_data_2)
  }
}

forest_plot = ggplot(All_data) + theme_bw() + 
  geom_segment(aes(x = UL, xend = LL, y = Group, yend = Group, 
                   color = Group), linewidth = 1) + 
  geom_point(aes(x = OR, y = Group, color = Group), size = 2) + 
  scale_color_manual(breaks = c("cdsrABP", "srABP", 
                                "cdABP", "Control"), 
                     values = c("#65ACEA", "#F9D057",
                                "#D05783", "grey25")) + 
  facet_wrap(~Variable, nrow = 3) + 
  theme(legend.position = 'top') + 
  geom_vline(xintercept = 1) + 
  xlab("Odds Ratio (95% Confidence Interval)") + 
  ylab(element_blank()) + 
  coord_cartesian(xlim = c(-0.75, 4.5)) + 
  geom_text(data = Case_Confounders_2, 
            aes(x = -0.5, y = Group, 
                label = paste("p=",scientific(P, 2), sep = "")),
            size = 3.5);forest_plot
jpeg("./Images/IBS_4/forest_plot.jpeg",
     res = 400, height = 2000, width = 2000)
forest_plot
dev.off()
