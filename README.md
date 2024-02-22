### Step 1: Phyloseq_Generation_IBS.R 
**Inputs:** 
1. Reported data from ABIS surveys and Swedish National Registry in .dta form 
2. Metadata for ABIS Infants 
3. Phyloseq object for ABIS Infant Stool samples 
4. CSV file for Age of diagnosis 

**Outputs:**
1. Sample Sheet of Related IBS data 
2. Phyloseq object for selected IBS subjects and controls 
3. Boxplot of age of diagnosis for selected IBS groups 

### Step 2: Chi-sq heterogeneity: Heterogeneous.R
**Inputs:** 
1. IBS Sample data
**Outputs:**
1. CSV of heterogeneous factors
2. Forest plot of odds Ratio

### Step 3: Beta diversity: BetaDiversity_Confounders.R
**Inputs:** 
1. IBS Phyloseq object 
**Outputs:**
2. CSV of ANOVA p values 

### Step 4: Alpha_Diversity.R
**Inputs:** 
1. IBS Phyloseq Object 
**Outputs:**
1. Boxplots of rarefied alpha diversity

### Step 5: Phyla Wilcoxon 
**Inputs:** 
**Outputs:**

### Step 6: PIME OOB: PIME_Prevalence.R
**Inputs:** 
1. IBS Phyloseq object
**Outputs:**
1. Boxplots of OOB from 50 iterations of Random Forest
2. CSV files of genera with MDA > 0

### Step 7: LEfSE and ALDEx2: Differential Abundance 
**Inputs:** 
1. IBS Phyloseq Object 
**Outputs:**
1. CSV File of LEFSE markers
2. CSV File of ALDEx2 markers
3. Dotplot of lefse and Aldex2 markers 

### Step 8: Pheatmap 
**Inputs:** 
1. IBS Phyloseq object
2. CSV of PIME markers
3. CSV of LEFSE genera
4. CSV of ALDEX genera 
**Outputs:**
1. Pheatmap of genera found in two or more analyses 

### Step 9: Picrust File Generation  
**Inputs:** 
1. IBS phyloseq object 
**Outputs:**
1. FASTA with identifiers (Genera_RowNumber) and ASV sequence
2. Count of total abundance per subject and identifier
    - picrust2 script:
        - conda activate picrust2
        - picrust2_pipeline.py -s IBS_TotalAbun.fasta -i IBS_count_TotalAbun.tsv -o picrust2_out_pipeline -p 2

### Step 10: Picrust Analysis: Picrust_Analysis.R
**Inputs:** 
1. IBS phyloseq object
2. pathways_out/path_abun_unstrat.tsv
3. pathways_out/metacyc_pathways_info.csv
**Outputs:**
1. Boxplot of relative abundance of predicted pathways
2. Boxplot of predicted pathways vs. heterogeneous factors 

### Step 11: Environmental Confounders 
**Inputs:** 
**Outputs:**

