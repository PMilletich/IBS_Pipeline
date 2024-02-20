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

### Step 2: Chi-sq heterogeneity 
**Inputs:** 
**Outputs:**

### Step 3: Beta diversity 
**Inputs:** 
**Outputs:**

### Step 4: Alpha Diversity 
**Inputs:** 
**Outputs:**

### Step 5: Phyla Wilcoxon 
**Inputs:** 
**Outputs:**

### Step 6: PIME OOB 
**Inputs:** 
**Outputs:**

### Step 7: LEfSE and ALDEx2 
**Inputs:** 
**Outputs:**

### Step 8: Pheatmap 
**Inputs:** 
**Outputs:**


### Step 9: Picrust File Generation  
**Inputs:** 
1. IBS phyloseq object 

**Outputs:**
1. FASTA with identifiers (Genera_RowNumber) and ASV sequence
2. Count of total abundance per subject and identifier
    - picrust2 script:
        - conda activate picrust2
        - picrust2_pipeline.py -s IBS_TotalAbun.fasta -i IBS_count_TotalAbun.tsv -o picrust2_out_pipeline -p 2

### Step 10: Picrust Analysis 
**Inputs:** 
1. IBS phyloseq object
2. KO_metagenome_out/pred_metagenome_unstrat.tsv
3. KO_metagenome_out/ko_info.csv
   - https://github.com/picrust/picrust2/tree/master/picrust2/default_files/description_mapfiles
     
**Outputs:**
1. csv of significant predictions
2. boxplots of significant predictions


### Step 11: Environmental Confounders 
**Inputs:** 
**Outputs:**

