library(haven) #https://cran.r-project.org/web/packages/haven/index.html
library(ggplot2) #https://ggplot2.tidyverse.org/
library(phyloseq) #https://bioconductor.org/packages/release/bioc/html/phyloseq.html
library(dplyr) #https://cran.r-project.org/web/packages/dplyr/index.html

setwd("~/Desktop/IBS/")

############################################################
#IBS Data management
############################################################
IBS_data = read_dta("Abdpain_microbiome.dta")
IBS_data$ID = IBS_data$abisnr
IBS_data <- data.frame(sapply(IBS_data, function(x) as.numeric(as.character(x))))

#Create compilation columns if there are any timepoint of abp
IBS_data_NA = IBS_data
IBS_data_NA[is.na(IBS_data_NA)] = 0
IBS_data_NA$srABP = apply(IBS_data_NA[,c("abp2","abp5","abp8","abp12","abp20")], 1, sum)
IBS_data_NA$srABP = ifelse(is.na(IBS_data_NA$srABP) == T, NA, 
                         ifelse(IBS_data_NA$srABP >= 1, 1, 0))

#Compile Functional Gastrointestinal diagnoses
IBS_data_NA$cdFG = apply(IBS_data_NA[,c("k58", "k590","k59","k591","k598","k599" )], 1, sum)
IBS_data_NA$cdFG = ifelse(is.na(IBS_data_NA$cdFG) == T, NA, 
                                  ifelse(IBS_data_NA$cdFG >= 1, 1, 0))

#Compile abp clinical diagnoses
IBS_data_NA$cdABP = apply(IBS_data_NA[,c("r10","r104" )], 1, sum)
IBS_data_NA$cdABP = ifelse(is.na(IBS_data_NA$cdABP) == T, NA, 
                                ifelse(IBS_data_NA$cdABP >= 1, 1, 0))

#Create column of total symptoms; focusing on abp
IBS_data_NA$IBS = rowSums(IBS_data_NA[,c("srABP","cdABP")], na.rm = T)
IBS_data_NA$IBS = ifelse(IBS_data_NA$IBS > 1, "cdsrMULT", 
                         ifelse(IBS_data_NA$srABP > 0, "srABP", 
                                ifelse(IBS_data_NA$cdABP > 0, "cdABP",
                                       ifelse(IBS_data_NA$cdFG > 0, "cdFG", "Control"))))


############################################################
#Merge with Sample data 
############################################################
Stool_data = read.csv("~/Desktop/ABIS/ABIS_1year_6_2023.csv")
PS = readRDS("~/Desktop/ABIS/RDP_2023.RDS")

Merge_data = merge(IBS_data_NA, Stool_data) #632
rownames(Merge_data) = Merge_data$ID
table(Merge_data$IBS)

#######
#Remove Celiac or Crohn's
table(Merge_data$Autoimmune_Disorder)
Merge_data = subset(Merge_data, grepl("Celiac", Merge_data$Autoimmune_Disorder) == F & 
                      grepl("Crohns", Merge_data$Autoimmune_Disorder) == F)
#617
rownames(Merge_data) = Merge_data$ID

##################################################################
#Remove Low counts; remove with less than 1000 raw reads
#Remove those missing qpcr data (copies_16s_per_gram_stool)
##################################################################
sample_data(PS) = sample_data(Merge_data)
samplesums = data.frame(sample_sums(PS))
samplesums$ID =rownames(samplesums)
summary(samplesums$sample_sums.PS.)
samplesums_1 = subset(samplesums, samplesums$sample_sums.PS. >= 1000)

Merge_data = subset(Merge_data,
                    is.na(Merge_data$copies_16s_per_gram_stool) == F &
                      Merge_data$ID %in% samplesums_1$ID)
sample_data(PS) = sample_data(Merge_data)
#[ 10583 taxa and 605 samples ]
##################################################################
#Age Subset 
##################################################################
Merge_data$Age.collected.month = round(Merge_data$Age.collected.month)  

ggplot(subset(Merge_data, Merge_data$IBS != "Control"),
       aes(x = Age.collected.month, fill = IBS)) + 
  geom_bar() + theme_bw()

Merge_data_Age = subset(Merge_data, Merge_data$Age.collected.month >= 10 & 
                          Merge_data$Age.collected.month <= 14)
#474
table(Merge_data_Age$IBS)
table(Merge_data_Age$Age.collected.month)

##################################################################
#Remove Functional cdFG
##################################################################
Merge_data_Age = subset(Merge_data_Age, Merge_data_Age$IBS != "cdFG")
#476
##################################################################
#Add Age of Diagnosis 
##################################################################
Diagnosis_age = read.csv("./CSV_Files/IBS_AgeDiag.csv")
sample_age = merge(Diagnosis_age, Merge_data_Age, all.y = T)

#' For those that have both cdABP and srABP,
#' Groups based on the youngest reported age 
#' If cdABP has no reported age, bin to srABP 
sample_age_mult = subset(sample_age, sample_age$IBS == "cdsrMULT")
cdABP_i = 0
srABP_i = 0
NA_i = 0
for (i in unique(sample_age_mult$ID)) {
  current_row = subset(sample_age_mult, sample_age_mult$ID == i)
  cd_age = current_row$Age_IBS
  sr_age = ifelse(current_row$abp2 == 1, 2, 
                ifelse(current_row$abp5 == 1, 5, 
                       ifelse(current_row$abp8 == 1, 8, 
                              ifelse(current_row$abp12 == 1, 12,
                                     ifelse(current_row$abp20== 1,18,NA)))))
  if (is.na(cd_age) == T & is.na(sr_age) == F) {
    sample_age[sample_age$ID == i, "Age_IBS"] = sr_age
    sample_age[sample_age$ID == i, "IBS"] = "srABP"
    srABP_i = srABP_i + 1
  } else if (is.na(cd_age) == T & is.na(sr_age) == T) {
    sample_age[sample_age$ID == i, "Age_IBS"] = NA
    sample_age[sample_age$ID == i, "IBS"] = "cdABP"
    NA_i = NA_i + 1
  } else if (cd_age < sr_age) {
    sample_age[sample_age$ID == i, "Age_IBS"] = cd_age
    sample_age[sample_age$ID == i, "IBS"] = "cdABP"
    cdABP_i = cdABP_i + 1 
  } else if (sr_age < cd_age) {
    sample_age[sample_age$ID == i, "Age_IBS"] = sr_age
    sample_age[sample_age$ID == i, "IBS"] = "srABP"
    srABP_i = srABP_i + 1
  } else {
    print(i)
  }
}
print(paste(cdABP_i, srABP_i, NA_i))
################################################################
#Create boxplot of age of diagnosis 
################################################################
sample_age_IBS = subset(sample_age, sample_age$IBS != "Control")
sample_age_IBS$IBS = "cdsrABP"
sample_age_IBS = rbind(sample_age_IBS, sample_age)

for (i in unique( sample_age_IBS$IBS)) {
  i_subset = subset(sample_age_IBS, sample_age_IBS$IBS == i)

    print(paste(i, nrow(i_subset),
                round(mean(i_subset$Age_IBS, na.rm = T),2), 
                round(sd(i_subset$Age_IBS, na.rm = T), 2)))

}

table(sample_age_IBS$IBS)

boxplot = ggplot(subset(sample_age_IBS, sample_age_IBS$IBS != "Control"), 
                 aes(x=round(Age_IBS), y = IBS, fill = IBS)) + 
  geom_boxplot() + theme_bw() + 
  xlab("Month of IBS Report") + 
  theme(legend.position = "top")+
  scale_fill_manual(breaks = c("cdsrABP", "srABP", "cdABP", "Control"), 
                     values = c("#65ACEA", "#F9D057", "#D05783", "grey25")); boxplot


rownames(sample_age) = sample_age$ID
sample_data(PS) = sample_data(sample_age)

############################################################
#Report data on ASVs and Genera 
############################################################
samplesums = data.frame(sample_sums(PS))
mean(samplesums[,1]); max(samplesums[,1])
ps.filtered = filter_taxa(PS, function(x) sum(x >= 1) > 1, TRUE)
ps.genus = tax_glom(ps.filtered, "Genus")

############################################################
#Save Files
############################################################
write.csv(sample_age, "IBS_1Year_IBS.csv", row.names = F)
saveRDS(ps.filtered, "IBS.RDS")

jpeg("./Images/IBS_4/Age_Diagnosis.jpeg", res = 400, 
     height = 2000, width = 2000)
boxplot
dev.off()
