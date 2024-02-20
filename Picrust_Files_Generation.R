#Libraries 
library(seqinr) #https://cran.r-project.org/web/packages/seqinr/index.html

#########################################################
#File input 
setwd("~/Desktop/IBS/")

OG_PS = readRDS("IBS.RDS")
Sample_data = data.frame(sample_data(OG_PS))
rownames(Sample_data) = Sample_data$ID

##########################################################
#Create unique identifiers for Genera ASVs
taxa_id = data.frame(tax_table(OG_PS))
taxa_id$ASV = rownames(taxa_id)
taxa_id = taxa_id[order(taxa_id$Genus),]
rownames(taxa_id) = 1:nrow(taxa_id)
taxa_id$Identifier = paste(taxa_id$Genus, rownames(taxa_id), sep = "_")
taxa_id = taxa_id[,c("ASV", "Identifier")]

###################
#Absolute vs Relative Abundance 
###################
ps.RA = transform_sample_counts(OG_PS, function(x) x / sum(x) )

qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
row.names(qpcr_data) = qpcr_data$ID
qpcr_data$ID = NULL

#Merge the Relative Abundance table and the qpcr table 
otu_RA = data.frame(otu_table(ps.RA))
otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
row.names(otu_RA_1) = otu_RA_1$Row.names
otu_RA_1$Row.names = NULL
i = 4
#Multiply all the Relabundance columns by the reads/g column
for(i in 2:length(colnames(otu_RA_1))) {
  otu_RA_1[is.na(otu_RA_1[,i]),i] = 0
  otu_RA_1[,i] <- suppressWarnings(as.integer(round(otu_RA_1[,1] * otu_RA_1[, i])))
}
ps.RA.qpcr_1 = ps.RA
otu_table(ps.RA.qpcr_1) = otu_table(otu_RA_1, taxa_are_rows = FALSE)

#########################################################
#create tsv table with new identifiers 
#########################################################
IBS_OTU = data.frame(t(otu_table(ps.RA.qpcr_1)))
IBS_OTU$ASV = rownames(IBS_OTU)
IBS_OTU = merge(IBS_OTU, taxa_id)
rownames(IBS_OTU) = IBS_OTU$Identifier
IBS_OTU$ASV = NULL;IBS_OTU$Identifier = NULL

write.table(IBS_OTU, "IBS_count_TotalAbun.tsv", sep = "\t")

############################################################
#Fastq File
############################################################
IBS_Tax = data.frame(tax_table(ps.RA.qpcr_1))
IBS_Tax$ASV = rownames(IBS_Tax)
IBS_Tax = merge(IBS_Tax, taxa_id)
rownames(IBS_Tax) = IBS_Tax$Identifier
#Ensure in same order as OTU Table
IBS_Tax = IBS_Tax[rownames(IBS_OTU),]
IBS_Tax = IBS_Tax[,c("Identifier", "ASV")] 
colnames(IBS_Tax) = c("name", "seq")
#create proper fastq headers
IBS_Tax$name = paste(">", IBS_Tax$name, sep = "")
IBS_Tax = subset(IBS_Tax, IBS_Tax$name != ">NA")

#Save as FASTA 
D <- do.call(rbind, lapply(seq(nrow(IBS_Tax)), function(i) t(IBS_Tax[i, ])))
write.table(D, file = "IBS_TotalAbun.fasta", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

