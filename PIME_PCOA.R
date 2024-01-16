#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(pime)
library(ggplot2)
library(ggpubr)

#################################
###File Import 
physeq_f = readRDS("./IBS.RDS")
physeq_f = tax_glom(physeq_f, "Genus")
sample = read.csv("./IBS_1Year.csv")
sample$abp_Combined = ifelse(is.na(sample$abp_Multiple) == T, NA, 
                             ifelse(sample$abp_Multiple == "No", "No", "Any"))
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
#Create dataframe to append top MDA genera for each four models 
All_MDA = data.frame()

##############################################
#Any IBS 
ps.RA.qpcr_Any = ps.RA.qpcr
set.seed(1)
pime.oob.error(ps.RA.qpcr_Any, "Any")
per_variable_obj_Any= pime.split.by.variable(ps.RA.qpcr_Any, "Any")
prevalences_Any=pime.prevalence(per_variable_obj_Any)
set.seed(42)
best.prev_Any=pime.best.prevalence(prevalences_Any, "Any")

#Most important Taxa
current_imp_Any=best.prev_Any$`Importance`$`Prevalence 65`
current_imp_Any$Group = "Any"
current_imp_Any = current_imp_Any[,c("Yes", "No", "MeanDecreaseAccuracy", 
                             "Genus", "Group")]
All_MDA = rbind(All_MDA, current_imp_Any)

#Prevalence Filtering
prevalence.current_Any= prevalences_Any$`65`
input_ord_Any = ordinate(prevalence.current_Any, "PCoA" , "binomial")

Any_PCOA= plot_ordination(prevalence.current_Any, input_ord_Any , 
                          color = "Any", 
                           shape = "Any")+
  stat_ellipse(aes(group = Any)) + 
  scale_color_manual(breaks = c("No", "Yes"), 
                     values = c("grey50", "red3")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")); Any_PCOA


##############################################
#ABP
abp_df = subset(sample, sample$abp_Combined == "Any" | 
                  sample$Any == "No")
table(sample$abp_Combined); table(abp_df$abp_Combined)

ps.RA.qpcr_abp = ps.RA.qpcr
sample_data(ps.RA.qpcr_abp) = sample_data(abp_df)
set.seed(1)
pime.oob.error(ps.RA.qpcr_abp, "abp_Combined")
per_variable_obj_abp= pime.split.by.variable(ps.RA.qpcr_abp, "abp_Combined")
prevalences_abp=pime.prevalence(per_variable_obj_abp)
set.seed(42)
best.prev_abp=pime.best.prevalence(prevalences_abp, "abp_Combined")

#Top Genera
current_imp_abp=best.prev_abp$`Importance`$`Prevalence 65`
current_imp_abp$Group = "abp"
current_imp_abp = current_imp_abp[,c("Any", "No", "MeanDecreaseAccuracy", 
                                     "Genus", "Group")]
colnames(current_imp_abp) = c("Yes", "No", "MeanDecreaseAccuracy", "Genus", "Group")
All_MDA = rbind(All_MDA, current_imp_abp)

prevalence.current_abp= prevalences_abp$`65`
input_ord_abp = ordinate(prevalence.current_abp, "PCoA" , "binomial")

abp_PCOA= plot_ordination(prevalence.current_abp, input_ord_abp , 
                          color = "abp_Combined", 
                          shape = "abp_Combined")+
  stat_ellipse(aes(group = abp_Combined)) + 
  scale_color_manual(breaks = c("Any", "No"),
                     values = c("purple3", "grey50")) +
  guides(color=guide_legend(title="Abdominal Pain\ntimepoints"), 
         shape=guide_legend(title="Abdominal Pain\ntimepoints")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")); abp_PCOA

##############################################
#k
K_Sample = subset(sample, sample$k_Combined == "Any" | 
                      sample$Any == "No")
ps.RA.qpcr_k = ps.RA.qpcr
sample_data(ps.RA.qpcr_k) = sample_data(K_Sample)
table(sample$k_Combined); table(K_Sample$k_Combined)

set.seed(1)
pime.oob.error(ps.RA.qpcr_k, "k_Combined")
per_variable_obj= pime.split.by.variable(ps.RA.qpcr_k, "k_Combined")
prevalences=pime.prevalence(per_variable_obj)
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "k_Combined")

current_imp=best.prev$`Importance`$`Prevalence 65`
current_imp$Group = "k"
current_imp = current_imp[,c("Any", "None", "MeanDecreaseAccuracy", "Genus", "Group")]
colnames(current_imp) = c("Yes", "No", "MeanDecreaseAccuracy", "Genus", "Group")
All_MDA = rbind(All_MDA, current_imp)

prevalence.current= prevalences$`65`
input_ord = ordinate(prevalence.current, "PCoA" , "binomial")

k_PCOA= plot_ordination(prevalence.current, input_ord , color = "k_Combined", 
                          shape = "k_Combined")+
  stat_ellipse(aes(group = k_Combined)) + 
  scale_color_manual(breaks = c("Any", "None"), 
                     values = c("forestgreen", "grey50")) + 
  guides(color=guide_legend(title="K Symptoms"), 
         shape=guide_legend(title="K Symptoms")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")); k_PCOA

##############################################
#r
r_Sample = subset(sample, sample$r_Combined == "Any" | 
                    sample$Any == "No")
ps.RA.qpcr_r = ps.RA.qpcr
sample_data(ps.RA.qpcr_r) = sample_data(r_Sample)
table(sample$r_Combined); table(r_Sample$r_Combined)

set.seed(1)
pime.oob.error(ps.RA.qpcr_r, "r_Combined")
per_variable_obj= pime.split.by.variable(ps.RA.qpcr_r, "r_Combined")
prevalences=pime.prevalence(per_variable_obj)
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "r_Combined")

current_imp=best.prev$`Importance`$`Prevalence 65`
current_imp$Group = "r"
current_imp = current_imp[,c("Any", "None", "MeanDecreaseAccuracy", "Genus", "Group")]
colnames(current_imp) = c("Yes", "No", "MeanDecreaseAccuracy", "Genus", "Group")
All_MDA = rbind(All_MDA, current_imp)

prevalence.current= prevalences$`65`
input_ord = ordinate(prevalence.current, "PCoA" , "binomial")

r_PCOA= plot_ordination(prevalence.current, input_ord , color = "r_Combined", 
                        shape = "r_Combined")+
  stat_ellipse(aes(group = r_Combined)) + 
  scale_color_manual(breaks = c("Any", "None"), 
                     values = c("blue3", "grey50")) + 
  guides(color=guide_legend(title="R Symptoms"), 
         shape=guide_legend(title="R Symptoms")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")); r_PCOA

jpeg("./Images/PCOA_4Groups_3.jpeg", res = 400, height = 4000, width = 6000)
ggarrange(Any_PCOA, 
          ggarrange(abp_PCOA, r_PCOA, k_PCOA, ncol = 3), 
          nrow = 2)
dev.off()

All_MDA_sign = subset(All_MDA, All_MDA$Yes > 0 & 
                        All_MDA$No > 0)

MDA_Count = data.frame(table(All_MDA_sign$Genus))
MDA_Count_1 = subset(MDA_Count, MDA_Count$Freq > 1)

write.csv(All_MDA_sign, "./CSV_Files/PIME_All.csv")
