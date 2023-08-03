#RNA-Seq CemiTool and network analysis pipeline 
#Valdes et. al 2023 Molecular Psychiatry Submission
#Load packages
library("edgeR")
library("limma")
library("data.table")
library("CEMiTool")
library("tmod")
library("fgsea")

#############################################################################
#CEMiTool Network Analysis Pipeline
#############################################################################
# load limma y object
load("~/PSEN1.PSEN2.APP.y.Rdata")
load("~/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")
samples <- y$samples
Psamples <- samples
Psamples$SampleName <- Psamples$sample
Psamples$Class <- Psamples$condition
sample_annot_expr <- Psamples[,c(11,12)]
rownames(sample_annot_expr) <- sample_annot_expr$SampleName
sample_annot_expr$Class <- paste('Cluster', sample_annot_expr$Class, sep='_')
sample_annot_expr$Class <- as.character(sample_annot_expr$Class)
sample_annot_expr$SampleName <- as.character(sample_annot_expr$SampleName)
counts_expr <- y$counts
rownames(counts_expr) <- y$genes$SYMBOL[y$genes$ENSEMBL %in% rownames(counts_expr)]
head(counts_expr)

#Create the cem object
cem <- cemitool(counts_expr, sample_annot_expr, 
                apply_vst = TRUE,filter=TRUE, filter_pval = 0.05,
                network_type = "signed")

#Get module hubs 
hubs <- get_hubs(cem, n=10, "adjacency")

#Get modules for only module genes in the modules
modules <- unique(cem@module[, 'modules'])
gene_sets <- lapply(modules, function(mod){
  return(cem@module[cem@module[, 'modules']==mod, 'genes'])
})
names(gene_sets) <- modules
cem.modules.gmt <- gene_sets

###################################
#2. Creation of Module 1 Networks
###################################

#Grab the M1 genes
M1_genes <- as.data.frame(gene_sets$M1)

#Rename column where names is "gene_sets$M1"
names(M1_genes)[names(M1_genes) == "gene_sets$M1"] <- "genesymbol1"

#Read in the interactions files 
base_dir <- '/Users/phoebevaldes/Desktop/CEMitool/Interactions'

#########################
#2a. EnCODE/ChEA for M1
#########################

# load ENCODE-ChEA consensus
ECC <- read.table(file.path(base_dir, "ENCODE_CHEA_CONSENSUS_FINAL.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Grab the unique TF terms
ECC_TF <- unique(ECC$gene1symbol)
ECC_TF_M1 <- c(ECC_TF,as.character(M1_genes$genesymbol1))
ECC_TF_M1_df <- as.data.frame(ECC_TF_M1)

#Rename the column
names(ECC_TF_M1_df)[names(ECC_TF_M1_df) == "ECC_TF_M1"] <- "genesymbol1"

##################
#2b. ReMap for M1
##################

# Load ReMap database
ReMap <- read.table(file.path(base_dir, "ReMap_FINAL.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Grab the unique TF terms
ReMap_TF <- unique(ReMap$gene1symbol)
ReMap_TF_M1 <- c(ReMap_TF,as.character(M1_genes$genesymbol1))
ReMap_TF_M1_df <- as.data.frame(ReMap_TF_M1)

#Rename the column
names(ReMap_TF_M1_df)[names(ReMap_TF_M1_df) == "ReMap_TF_M1"] <- "genesymbol1"

#######################################################
#2c. Get PPI interactions from stringDB for Module 1
#######################################################

#Read in string DB object
string_db <- read.table(file.path(base_dir, "stringdb_exdb.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Find dimensions of PPI edges object
string_db_dim <- dim(string_db)

#Create new stringDB object with color code of '1'
new_string_db <- string_db %>% 
  mutate(Color_Code = '1')

#Transpose the unfiltered expression matrix
new_string_db <- as.data.frame(t(as.matrix(new_string_db)))
View(new_string_db)

#Rename column names in numeric matrix
#Rename the columns
(setattr(new_string_db, "row.names", c("gene1symbol", "gene2symbol", "Color_Code")))
new_string_db <- t(as.data.frame(new_string_db))

new_string_db_unique <- as.data.frame(unique(new_string_db))

#Subset for only M1 genes for PPI interactions + neighboring genes
new_stringdb_M1 <- new_string_db_unique[new_string_db_unique$gene1symbol %in% M1_genes$genesymbol1,]

##################################################################
#2d. Get TF-gene interactions from ENCODE/ChEA and ReMap for Module 1
##################################################################

#Subset for only M1 genes for TF-gene interactions 
new_ECC_M1 <- ECC[ECC$gene1symbol %in% M1_genes$genesymbol1,]
new_ECC_M1 <- new_ECC_M1[new_ECC_M1$gene2symbol %in% M1_genes$genesymbol1,]

#Create new ECC subset object with color code of '2'
new_ECC_M1 <- new_ECC_M1 %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ECC_M1 <- unique(new_ECC_M1)

#Create new dataframe with one table on top of the other
new_string_db_ECC_M1 <- rbind(new_stringdb_M1, new_ECC_M1)

#Subset Module 1 TF's from ECC interactions
ECC_subset <- 
  subset(ECC, gene1symbol == c("FOXM1", "POU5F1", "TRIM28", "TP63"))

#Create new ECC subset object with color code of '2'
new_ECC_subset <- ECC_subset %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ECC_subset <- unique(new_ECC_subset)

#Get ECC subset with only M1 genes for genesymbol2
test_ECC <- new_ECC_subset[new_ECC_subset$gene2symbol %in% M1_genes$genesymbol1,]


#Create another new dataframe with another table on top of the other
new_string_db_ECC_M1 <- rbind(new_string_db_ECC_M1, new_ECC_subset)

#Create new ECC subset object with color code of '2'
new_ECC_M1 <- new_ECC_M1 %>% 
  mutate(Color_Code = '2')

#Create another new dataframe with another table on top of the other
new_ECC_M1_df <- rbind(new_ECC_M1, test_ECC)
new_ECC_M1_df <- unique(new_ECC_M1_df)

#Subset for only M1 genes for ReMap TF-gene interactions 
new_ReMap_M1 <- ReMap[ReMap$gene1symbol %in% M1_genes$genesymbol1,]
new_ReMap_M1 <- new_ReMap_M1[new_ReMap_M1$gene2symbol %in% M1_genes$genesymbol1,]

#Create new ECC subset object with color code of '2'
new_ReMap_M1 <- new_ReMap_M1 %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ReMap_M1 <- unique(new_ReMap_M1)

#Subset Module 1 TF's from ReMap interactions
ReMap_subset <- 
  subset(ReMap, gene1symbol == c("PCGF2"))

#Create new ECC subset object with color code of '2'
new_ReMap_subset <- ReMap_subset %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ReMap_subset <- unique(new_ReMap_subset)

#Get ReMap subset with only M1 genes for genesymbol2
test_ReMap <- new_ReMap_subset[new_ReMap_subset$gene2symbol %in% M1_genes$genesymbol1,]

#Create new dataframe with one table on top of the other
new_string_db_ReMap <- rbind(new_string_db_unique, new_ReMap_subset)

#Produce test objects that worked for stringDB interactions (find M1 genes in ReMap + stringDB interactions)
test_string_db_ReMap <- new_string_db_ReMap[new_string_db_ReMap$gene1symbol %in% M1_genes$genesymbol1,]
View(test_string_db_ReMap)

#Create another new dataframe with another table on top of the other
new_ReMap_M1_df <- rbind(new_ReMap_M1, test_ReMap)
new_ReMap_M1_df <- unique(new_ReMap_M1_df)

#Create new dataframe with ECC and ReMap interactions
new_ECC_ReMap_M1 <- rbind(new_ECC_M1_df, new_ReMap_M1_df)

################################################################
#2e. Enrichment for PPI interactions for Module 1 (APP vs. NDC)
################################################################

#Included M1 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_V717I <- cbind(gene2symbol = rownames(reslimma_V717I), reslimma_V717I)
rownames(reslimma_V717I) <- 1:nrow(reslimma_V717I)

#Subset M1 genes into only stringDB 
#Merge two data frames by ID
M1_string_db_ECC_merge_V717I_df <- merge(new_stringdb_M1,reslimma_V717I,by="gene2symbol")

#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M1_string_db_ECC_merge_V717I_df, file="~/M1/M1_string_db_merge_V717I.csv", row.names = F) 

#Note: checking which gene1symbols are M1 genes
test_genesymbol1 <- intersect(M1_string_db_ECC_merge_V717I_df$gene1symbol, M1_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M1/M1_string_db_genesymbol1_intersect_PPI_V717I.csv", row.names = F) 

#Note: checking which gene2symbols are M1 genes
test_genesymbol2 <- intersect(M1_string_db_ECC_merge_V717I_df$gene2symbol, M1_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M1/M1_string_db_genesymbol2_intersect_PPI_V717I.csv", row.names = F)

##############################################################################################
#2f. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 1 (APP vs. NDC)
##############################################################################################

#Included M1 genes in ECC + ReMap interactions + additional TF's with only
#M1 genes as target genes

#Subset M1 genes into only stringDB 
#Merge two data frames by ID
M1_string_db_ECC_merge_V717I_df <- merge(new_ECC_ReMap_M1,reslimma_V717I,by="gene2symbol")

#Note: checking which gene1symbols are M1 genes
test_genesymbol1 <- intersect(M1_string_db_ECC_merge_V717I_df$gene1symbol, M1_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M1/M1_string_db_genesymbol1_intersect_ECC_ReMap_V717I.csv", row.names = F) 

#Note: checking which gene2symbols are M1 genes
test_genesymbol2 <- intersect(M1_string_db_ECC_merge_V717I_df$gene2symbol, M1_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M1/M1_string_db_genesymbol2_intersect_ECC_ReMap_V717I.csv", row.names = F) 

#New Table with EDGES information - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M1_string_db_ECC_merge_V717I_df, file="~/M1/M1_ECC_ReMap_merge_V717I.csv", row.names = F) 

####################################################################
#2g. Enrichment for PPI interactions for Module 1 (PSEN1 vs. NDC)
####################################################################

#Included M1 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_A79V <- cbind(gene2symbol = rownames(reslimma_A79V), reslimma_A79V)
rownames(reslimma_A79V) <- 1:nrow(reslimma_A79V)

#Subset M1 genes into only stringDB 
#Merge two data frames by ID
M1_string_db_ECC_merge_A79V_df <- merge(new_stringdb_M1,reslimma_A79V,by="gene2symbol")

#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M1_string_db_ECC_merge_A79V_df, file="~/M1/M1_string_db_merge_A79V.csv", row.names = F) 

#Note: checking which gene1symbols are M1 genes
test_genesymbol1_A79V <- intersect(M1_string_db_ECC_merge_A79V_df$gene1symbol, M1_genes$genesymbol1)  
View(test_genesymbol1_A79V)

#Convert character to data table 
test_genesymbol1_A79V <- as.data.table(test_genesymbol1_A79V)
write.csv(test_genesymbol1_A79V, file="~/M1/M1_string_db_genesymbol1_intersect_PPI_A79V.csv", row.names = F) 

#Note: checking which gene2symbols are M1 genes
test_genesymbol2_A79V <- intersect(M1_string_db_ECC_merge_A79V_df$gene2symbol, M1_genes$genesymbol1)  
View(test_genesymbol2_A79V)

#Convert character to data table 
test_genesymbol2_A79V <- as.data.table(test_genesymbol2_A79V)
write.csv(test_genesymbol2_A79V, file="~/M1/M1_string_db_genesymbol2_intersect_PPI_A79V.csv", row.names = F) 

###############################################################################################
#2h. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 1 (PSEN1 vs. NDC)
###############################################################################################

#Included M1 genes in ECC + ReMap interactions + additional TF's with only
#M1 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_A79V <- cbind(gene2symbol = rownames(reslimma_A79V), reslimma_A79V)
rownames(reslimma_A79V) <- 1:nrow(reslimma_A79V)

#Subset M1 genes into only stringDB 
#Merge two data frames by ID
M1_string_db_ECC_merge_A79V_df <- merge(new_ECC_ReMap_M1,reslimma_A79V,by="gene2symbol")


#Note: checking which gene1symbols are M1 genes
test_genesymbol1 <- intersect(M1_string_db_ECC_merge_A79V_df$gene1symbol, M1_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M1/M1_string_db_genesymbol1_intersect_ECC_ReMap_A79V.csv", row.names = F) 

#Note: checking which gene2symbols are M1 genes
test_genesymbol2 <- intersect(M1_string_db_ECC_merge_A79V_df$gene2symbol, M1_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M1/M1_string_db_genesymbol2_intersect_ECC_ReMap_A79V.csv", row.names = F) 

#New Table with EDGES information - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M1_string_db_ECC_merge_A79V_df, file="~/M1/M1_ECC_ReMap_merge_A79V.csv", row.names = F) 

#################################################################
#2i. Enrichment for PPI interactions for Module 1 (PSEN2 vs. NDC)
#################################################################

#Included M1 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_N141I <- cbind(gene2symbol = rownames(reslimma_N141I), reslimma_N141I)
rownames(reslimma_N141I) <- 1:nrow(reslimma_N141I)

#Subset M1 genes into only stringDB 
#Merge two data frames by ID
M1_string_db_ECC_merge_N141I_df <- merge(new_stringdb_M1,reslimma_N141I,by="gene2symbol")

#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M1_string_db_ECC_merge_N141I_df, file="~/M1/M1_string_db_merge_N141I_v3.csv", row.names = F) 

#Note: checking which gene1symbols are M1 genes
test_genesymbol1_N141I <- intersect(M1_string_db_ECC_merge_N141I_df$gene1symbol, M1_genes$genesymbol1)  
View(test_genesymbol1_N141I)

#Convert character to data table 
test_genesymbol1_N141I <- as.data.table(test_genesymbol1_N141I)
write.csv(test_genesymbol1_N141I, file="~/M1/M1_string_db_genesymbol1_intersect_PPI_N141I.csv", row.names = F) 

#Note: checking which gene2symbols are M1 genes
test_genesymbol2_N141I <- intersect(M1_string_db_ECC_merge_N141I_df$gene2symbol, M1_genes$genesymbol1)  
View(test_genesymbol2_N141I)

#Convert character to data table 
test_genesymbol2_N141I <- as.data.table(test_genesymbol2_N141I)
write.csv(test_genesymbol2_N141I, file="~/M1/M1_string_db_genesymbol2_intersect_PPI_N141I.csv", row.names = F) 

###############################################################################################
#2j. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 1 (PSEN2 vs. NDC)
###############################################################################################

#Included M1 genes in ECC + ReMap interactions + additional TF's with only
#M1 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_N141I <- cbind(gene2symbol = rownames(reslimma_N141I), reslimma_N141I)
rownames(reslimma_N141I) <- 1:nrow(reslimma_N141I)

#Subset M1 genes into only stringDB 
#Merge two data frames by ID
M1_string_db_ECC_merge_N141I_df <- merge(new_ECC_ReMap_M1,reslimma_N141I,by="gene2symbol")


#Note: checking which gene1symbols are M1 genes
test_genesymbol1 <- intersect(M1_string_db_ECC_merge_N141I_df$gene1symbol, M1_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M1/M1_string_db_genesymbol1_intersect_ECC_ReMap_N141I.csv", row.names = F) 

#Note: checking which gene2symbols are M1 genes
test_genesymbol2 <- intersect(M1_string_db_ECC_merge_N141I_df$gene2symbol, M1_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M1/M1_string_db_genesymbol2_intersect_ECC_ReMap_N141I.csv", row.names = F) 

#New Table with EDGES information  - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M1_string_db_ECC_merge_N141I_df, file="~/M1/M1_ECC_ReMap_merge_N141I.csv", row.names = F) 


#####################################################################################

###################################
#3. Creation of Module 3 Networks
###################################

#Grab the M3 genes
M3_genes <- as.data.frame(gene_sets$M3)

write.csv(M3_genes, file="~/M3/M3_genes.csv", row.names = F) 


#Rename column where names is "gene_sets$M3"
#Source: https://www.datanovia.com/en/lessons/rename-data-frame-columns-in-r/
names(M3_genes)[names(M3_genes) == "gene_sets$M3"] <- "genesymbol1"

########################
#3a. EnCODE/ChEA for M3
########################

ECC <- read.table(file.path(base_dir, "ENCODE_CHEA_CONSENSUS_FINAL.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Grab the unique TF terms
ECC_TF <- unique(ECC$gene1symbol)
ECC_TF_M3 <- c(ECC_TF,as.character(M3_genes$genesymbol1))
ECC_TF_M3_df <- as.data.frame(ECC_TF_M3)

#Rename the column
names(ECC_TF_M3_df)[names(ECC_TF_M3_df) == "ECC_TF_M3"] <- "genesymbol1"

##################
#3b. ReMap for M3
##################

ReMap <- read.table(file.path(base_dir, "ReMap_FINAL.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Grab the unique TF terms
ReMap_TF <- unique(ReMap$gene1symbol)
ReMap_TF_M3 <- c(ReMap_TF,as.character(M3_genes$genesymbol1))
ReMap_TF_M3_df <- as.data.frame(ReMap_TF_M3)

#Rename the column
names(ReMap_TF_M3_df)[names(ReMap_TF_M3_df) == "ReMap_TF_M3"] <- "genesymbol1"

#########################################################################################################

#Read in string DB object
string_db <- read.table(file.path(base_dir, "stringdb_exdb.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Find dimensions of PPI edges object
string_db_dim <- dim(string_db)

#Create new stringDB object with color code of '1'
new_string_db <- string_db %>% 
  mutate(Color_Code = '1')

#Transpose the unfiltered expression matrix
new_string_db <- as.data.frame(t(as.matrix(new_string_db)))
View(new_string_db)

#Rename column names in numeric matrix
#Rename the columns
(setattr(new_string_db, "row.names", c("gene1symbol", "gene2symbol", "Color_Code")))
new_string_db <- t(as.data.frame(new_string_db))

new_string_db_unique <- as.data.frame(unique(new_string_db))

#Subset for only M3 genes for PPI interactions + neighboring genes
new_stringdb_M3 <- new_string_db_unique[new_string_db_unique$gene1symbol %in% M3_genes$genesymbol1,]


#Subset for only M3 genes for TF-gene interactions 
new_ECC_M3 <- ECC[ECC$gene1symbol %in% M3_genes$genesymbol1,]
new_ECC_M3 <- new_ECC_M3[new_ECC_M3$gene2symbol %in% M3_genes$genesymbol1,]

#Create new ECC subset object with color code of '2'
new_ECC_M3 <- new_ECC_M3 %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ECC_M3 <- unique(new_ECC_M3)

#Create new dataframe with one table on top of the other
new_string_db_ECC_M3 <- rbind(new_stringdb_M3, new_ECC_M3)

#Subset Module 3 TF's from ECC interactions
ECC_subset <- 
  subset(ECC, gene1symbol == c("EZH2", "SALL4", "MYOD1", "ESR1", "SUZ12", "TRIM28"))

#Create new ECC subset object with color code of '2'
new_ECC_subset <- ECC_subset %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ECC_subset <- unique(new_ECC_subset)

#Get ECC subset with only M3 genes for genesymbol2
test_ECC <- new_ECC_subset[new_ECC_subset$gene2symbol %in% M3_genes$genesymbol1,]

#Create another new dataframe with another table on top of the other
new_string_db_ECC_M3 <- rbind(new_string_db_ECC_M3, new_ECC_subset)


#Create new ECC subset object with color code of '2'
new_ECC_M3 <- new_ECC_M3 %>% 
  mutate(Color_Code = '2')

#Create another new dataframe with another table on top of the other
new_ECC_M3_df <- rbind(new_ECC_M3, test_ECC)
new_ECC_M3_df <- unique(new_ECC_M3_df)

#Subset for only M3 genes for ReMap TF-gene interactions 
new_ReMap_M3 <- ReMap[ReMap$gene1symbol %in% M3_genes$genesymbol1,]
new_ReMap_M3 <- new_ReMap_M3[new_ReMap_M3$gene2symbol %in% M3_genes$genesymbol1,]

#Create new ECC subset object with color code of '2'
new_ReMap_M3 <- new_ReMap_M3 %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ReMap_M3 <- unique(new_ReMap_M3)

#Subset Module 3 TF's from ReMap interactions
ReMap_subset <- 
  subset(ReMap, gene1symbol == c("PCGF2", "NANOG", "TBXT", "SIX2"))

#Create new ECC subset object with color code of '2'
new_ReMap_subset <- ReMap_subset %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ReMap_subset <- unique(new_ReMap_subset)

#Get ReMap subset with only M3 genes for genesymbol2
test_ReMap <- new_ReMap_subset[new_ReMap_subset$gene2symbol %in% M3_genes$genesymbol1,]

#Create new dataframe with one table on top of the other
new_string_db_ReMap <- rbind(new_string_db_unique, new_ReMap_subset)

#Produce test objects that worked for stringDB interactions (find M3 genes in ReMap + stringDB interactions)
test_string_db_ReMap <- new_string_db_ReMap[new_string_db_ReMap$gene1symbol %in% M3_genes$genesymbol1,]
View(test_string_db_ReMap)


#Create another new dataframe with another table on top of the other
new_ReMap_M3_df <- rbind(new_ReMap_M3, test_ReMap)
new_ReMap_M3_df <- unique(new_ReMap_M3_df)

#Create new dataframe with ECC and ReMap interactions
new_ECC_ReMap_M3 <- rbind(new_ECC_M3_df, new_ReMap_M3_df)

##################################################################
#3c. Enrichment for PPI interactions for Module 3 (APP vs. NDC)
##################################################################

#Included M3 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_V717I <- cbind(gene2symbol = rownames(reslimma_V717I), reslimma_V717I)
rownames(reslimma_V717I) <- 1:nrow(reslimma_V717I)

#Subset M3 genes into only stringDB 
#Merge two data frames by ID
M3_string_db_ECC_merge_V717I_df <- merge(new_stringdb_M3,reslimma_V717I,by="gene2symbol")


#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M3_string_db_ECC_merge_V717I_df, file="~/M3/M3_string_db_merge_V717I_v3.csv", row.names = F) 

#Note: checking which gene1symbols are M3 genes
test_genesymbol1 <- intersect(M3_string_db_ECC_merge_V717I_df$gene1symbol, M3_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M3/M3_string_db_genesymbol1_intersect_PPI_V717I.csv", row.names = F) 

#Note: checking which gene2symbols are M3 genes
test_genesymbol2 <- intersect(M3_string_db_ECC_merge_V717I_df$gene2symbol, M3_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M3/M3_string_db_genesymbol2_intersect_PPI_V717I.csv", row.names = F)


#############################################################################################
#3d. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 3 (APP vs. NDC)
##############################################################################################

#Included M3 genes in ECC + ReMap interactions + additional TF's with only
#M3 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_V717I <- cbind(gene2symbol = rownames(reslimma_V717I), reslimma_V717I)
rownames(reslimma_V717I) <- 1:nrow(reslimma_V717I)

#Subset M3 genes into only stringDB 
#Merge two data frames by ID
M3_string_db_ECC_merge_V717I_df <- merge(new_ECC_ReMap_M3,reslimma_V717I,by="gene2symbol")


#Note: checking which gene1symbols are M3 genes
test_genesymbol1 <- intersect(M3_string_db_ECC_merge_V717I_df$gene1symbol, M3_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M3/M3_string_db_genesymbol1_intersect_ECC_ReMap_V717I.csv", row.names = F) 

#Note: checking which gene2symbols are M3 genes
test_genesymbol2 <- intersect(M3_string_db_ECC_merge_V717I_df$gene2symbol, M3_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M3/M3_string_db_genesymbol2_intersect_ECC_ReMap_V717I.csv", row.names = F) 

#New Table with EDGES information - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M3_string_db_ECC_merge_V717I_df, file="~/M3/M3_ECC_ReMap_merge_V717I.csv", row.names = F) 


#################################################################
#3f. Enrichment for PPI interactions for Module 3 (PSEN1 vs. NDC)
#################################################################

#Included M3 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_A79V <- cbind(gene2symbol = rownames(reslimma_A79V), reslimma_A79V)
rownames(reslimma_A79V) <- 1:nrow(reslimma_A79V)

#Subset M3 genes into only stringDB 
#Merge two data frames by ID
M3_string_db_ECC_merge_A79V_df <- merge(new_stringdb_M3,reslimma_A79V,by="gene2symbol")


#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M3_string_db_ECC_merge_A79V_df, file="~/M3/M3_string_db_merge_A79V_v3.csv", row.names = F) 

#Note: checking which gene1symbols are M3 genes
test_genesymbol1_A79V <- intersect(M3_string_db_ECC_merge_A79V_df$gene1symbol, M3_genes$genesymbol1)  
View(test_genesymbol1_A79V)

#Convert character to data table 
test_genesymbol1_A79V <- as.data.table(test_genesymbol1_A79V)
write.csv(test_genesymbol1_A79V, file="~/M3/M3_string_db_genesymbol1_intersect_PPI_A79V.csv", row.names = F) 

#Note: checking which gene2symbols are M3 genes
test_genesymbol2_A79V <- intersect(M3_string_db_ECC_merge_A79V_df$gene2symbol, M3_genes$genesymbol1)  
View(test_genesymbol2_A79V)

#Convert character to data table
test_genesymbol2_A79V <- as.data.table(test_genesymbol2_A79V)
write.csv(test_genesymbol2_A79V, file="~/M3/M3_string_db_genesymbol2_intersect_PPI_A79V.csv", row.names = F) 

###############################################################################################
#3g. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 3 (PSEN1 vs. NDC)
###############################################################################################

#Included M3 genes in ECC + ReMap interactions + additional TF's with only
#M3 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_A79V <- cbind(gene2symbol = rownames(reslimma_A79V), reslimma_A79V)
rownames(reslimma_A79V) <- 1:nrow(reslimma_A79V)

#Subset M3 genes into only stringDB 
#Merge two data frames by ID
M3_string_db_ECC_merge_A79V_df <- merge(new_ECC_ReMap_M3,reslimma_A79V,by="gene2symbol")


#Note: checking which gene1symbols are M3 genes
test_genesymbol1 <- intersect(M3_string_db_ECC_merge_A79V_df$gene1symbol, M3_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M3/M3_string_db_genesymbol1_intersect_ECC_ReMap_A79V.csv", row.names = F) 

#Note: checking which gene2symbols are M3 genes
test_genesymbol2 <- intersect(M3_string_db_ECC_merge_A79V_df$gene2symbol, M3_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M3/M3_string_db_genesymbol2_intersect_ECC_ReMap_A79V.csv", row.names = F) 

#New Table with EDGES information - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M3_string_db_ECC_merge_A79V_df, file="~/M3/M3_ECC_ReMap_merge_A79V.csv", row.names = F) 


##################################################################
#3h. Enrichment for PPI interactions for Module 3 (PSEN2 vs. NDC)
##################################################################

#Included M3 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_N141I <- cbind(gene2symbol = rownames(reslimma_N141I), reslimma_N141I)
rownames(reslimma_N141I) <- 1:nrow(reslimma_N141I)

#Subset M3 genes into only stringDB 
#Merge two data frames by ID
M3_string_db_ECC_merge_N141I_df <- merge(new_stringdb_M3,reslimma_N141I,by="gene2symbol")

#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M3_string_db_ECC_merge_N141I_df, file="~/M3/M3_string_db_merge_N141I_v3.csv", row.names = F) 

#Note: checking which gene1symbols are M3 genes
test_genesymbol1_N141I <- intersect(M3_string_db_ECC_merge_N141I_df$gene1symbol, M3_genes$genesymbol1)  
View(test_genesymbol1_N141I)

#Convert character to data table 
test_genesymbol1_N141I <- as.data.table(test_genesymbol1_N141I)
write.csv(test_genesymbol1_N141I, file="~/M3/M3_string_db_genesymbol1_intersect_PPI_N141I.csv", row.names = F) 

#Note: checking which gene2symbols are M3 genes
test_genesymbol2_N141I <- intersect(M3_string_db_ECC_merge_N141I_df$gene2symbol, M3_genes$genesymbol1)  
View(test_genesymbol2_N141I)

#Convert character to data table 
test_genesymbol2_N141I <- as.data.table(test_genesymbol2_N141I)
write.csv(test_genesymbol2_N141I, file="~/M3/M3_string_db_genesymbol2_intersect_PPI_N141I.csv", row.names = F) 

##############################################################################################
#3i. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 3 (PSEN2 vs. NDC)
###############################################################################################

#Included M3 genes in ECC + ReMap interactions + additional TF's with only
#M3 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_N141I <- cbind(gene2symbol = rownames(reslimma_N141I), reslimma_N141I)
rownames(reslimma_N141I) <- 1:nrow(reslimma_N141I)

#Subset M3 genes into only stringDB 
#Merge two data frames by ID
M3_string_db_ECC_merge_N141I_df <- merge(new_ECC_ReMap_M3,reslimma_N141I,by="gene2symbol")


#Note: checking which gene1symbols are M3 genes
test_genesymbol1 <- intersect(M3_string_db_ECC_merge_N141I_df$gene1symbol, M3_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M3/M3_string_db_genesymbol1_intersect_ECC_ReMap_N141I.csv", row.names = F) 

#Note: checking which gene2symbols are M3 genes
test_genesymbol2 <- intersect(M3_string_db_ECC_merge_N141I_df$gene2symbol, M3_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M3/M3_string_db_genesymbol2_intersect_ECC_ReMap_N141I.csv", row.names = F) 

#New Table with EDGES information  - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M3_string_db_ECC_merge_N141I_df, file="~/M3/M3_ECC_ReMap_merge_N141I.csv", row.names = F) 

#######################################################################################

##################################
#4. Creation of Module 4 Networks
##################################

#Grab the M4 genes
M4_genes <- as.data.frame(gene_sets$M4)

write.csv(M4_genes, file="~/M4/M4_genes.csv", row.names = F) 


#Rename column where names is "gene_sets$M3"
#Source: https://www.datanovia.com/en/lessons/rename-data-frame-columns-in-r/
names(M4_genes)[names(M4_genes) == "gene_sets$M4"] <- "genesymbol1"

########################
#4a. EnCODE/ChEA for M4
########################

ECC <- read.table(file.path(base_dir, "ENCODE_CHEA_CONSENSUS_FINAL.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Grab the unique TF terms
ECC_TF <- unique(ECC$gene1symbol)
ECC_TF_M4 <- c(ECC_TF,as.character(M4_genes$genesymbol1))
ECC_TF_M4_df <- as.data.frame(ECC_TF_M4)

#Rename the column
names(ECC_TF_M4_df)[names(ECC_TF_M4_df) == "ECC_TF_M4"] <- "genesymbol1"

##################
#4b. ReMap for M4
##################

ReMap <- read.table(file.path(base_dir, "ReMap_FINAL.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Grab the unique TF terms
ReMap_TF <- unique(ReMap$gene1symbol)
ReMap_TF_M4 <- c(ReMap_TF,as.character(M4_genes$genesymbol1))
ReMap_TF_M4_df <- as.data.frame(ReMap_TF_M4)

#Rename the column
names(ReMap_TF_M4_df)[names(ReMap_TF_M4_df) == "ReMap_TF_M4"] <- "genesymbol1"

#########################################################################################################

#Read in string DB object
string_db <- read.table(file.path(base_dir, "stringdb_exdb.tsv"), header = TRUE, stringsAsFactors=FALSE)

#Find dimensions of PPI edges object
string_db_dim <- dim(string_db)

#Create new stringDB object with color code of '1'
new_string_db <- string_db %>% 
  mutate(Color_Code = '1')

#Transpose the unfiltered expression matrix
new_string_db <- as.data.frame(t(as.matrix(new_string_db)))
View(new_string_db)

#Rename column names in numeric matrix
#Rename the columns
(setattr(new_string_db, "row.names", c("gene1symbol", "gene2symbol", "Color_Code")))
new_string_db <- t(as.data.frame(new_string_db))

new_string_db_unique <- as.data.frame(unique(new_string_db))

#Subset for only M4 genes for PPI interactions + neighboring genes
new_stringdb_M4 <- new_string_db_unique[new_string_db_unique$gene1symbol %in% M4_genes$genesymbol1,]

#Subset for only M4 genes for TF-gene interactions 
new_ECC_M4 <- ECC[ECC$gene1symbol %in% M4_genes$genesymbol1,]
new_ECC_M4 <- new_ECC_M4[new_ECC_M4$gene2symbol %in% M4_genes$genesymbol1,]

#Create new ECC subset object with color code of '2'
new_ECC_M4 <- new_ECC_M4 %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ECC_M4 <- unique(new_ECC_M4)

#Create new dataframe with one table on top of the other
new_string_db_ECC_M4 <- rbind(new_stringdb_M4, new_ECC_M4)

#Subset Module 4 TF's from ECC interactions
ECC_subset <- 
  subset(ECC, gene1symbol == c("SUZ12", "TRIM28"))

#Create new ECC subset object with color code of '2'
new_ECC_subset <- ECC_subset %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ECC_subset <- unique(new_ECC_subset)

#Get ECC subset with only M4 genes for genesymbol2
test_ECC <- new_ECC_subset[new_ECC_subset$gene2symbol %in% M4_genes$genesymbol1,]

#Create another new dataframe with another table on top of the other
new_string_db_ECC_M4 <- rbind(new_string_db_ECC_M4, new_ECC_subset)

#Create new ECC subset object with color code of '2'
new_ECC_M4 <- new_ECC_M4 %>% 
  mutate(Color_Code = '2')

#Create another new dataframe with another table on top of the other
new_ECC_M4_df <- rbind(new_ECC_M4, test_ECC)
new_ECC_M4_df <- unique(new_ECC_M4_df)

#Subset for only M4 genes for ReMap TF-gene interactions 
new_ReMap_M4 <- ReMap[ReMap$gene1symbol %in% M4_genes$genesymbol1,]
new_ReMap_M4 <- new_ReMap_M4[new_ReMap_M4$gene2symbol %in% M4_genes$genesymbol1,]

#Create new ECC subset object with color code of '2'
new_ReMap_M4 <- new_ReMap_M4 %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ReMap_M4 <- unique(new_ReMap_M4)

#Subset Module 4 TF's from ReMap interactions
ReMap_subset <- 
  subset(ReMap, gene1symbol == c("PCGF2", "NANOG"))

#Create new ECC subset object with color code of '2'
new_ReMap_subset <- ReMap_subset %>% 
  mutate(Color_Code = '2')

#Get the unique terms
new_ReMap_subset <- unique(new_ReMap_subset)

#Get ReMap subset with only M4 genes for genesymbol2
test_ReMap <- new_ReMap_subset[new_ReMap_subset$gene2symbol %in% M4_genes$genesymbol1,]

#Create new dataframe with one table on top of the other
new_string_db_ReMap <- rbind(new_string_db_unique, new_ReMap_subset)

#Produce test objects that worked for stringDB interactions (find M4 genes in ReMap + stringDB interactions)
test_string_db_ReMap <- new_string_db_ReMap[new_string_db_ReMap$gene1symbol %in% M4_genes$genesymbol1,]
View(test_string_db_ReMap)

#Create another new dataframe with another table on top of the other
new_ReMap_M4_df <- rbind(new_ReMap_M4, test_ReMap)
new_ReMap_M4_df <- unique(new_ReMap_M4_df)

#Create new dataframe with ECC and ReMap interactions
new_ECC_ReMap_M4 <- rbind(new_ECC_M4_df, new_ReMap_M4_df)

##################################################################
#4c. Enrichment for PPI interactions for Module 4 (APP vs. NDC)
##################################################################

#Included M4 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_V717I <- cbind(gene2symbol = rownames(reslimma_V717I), reslimma_V717I)
rownames(reslimma_V717I) <- 1:nrow(reslimma_V717I)

#Subset M4 genes into only stringDB 
#Merge two data frames by ID
M4_string_db_ECC_merge_V717I_df <- merge(new_stringdb_M4,reslimma_V717I,by="gene2symbol")

#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M4_string_db_ECC_merge_V717I_df, file="~/M4/M4_string_db_merge_V717I_v3.csv", row.names = F) 

#Note: checking which gene1symbols are M4 genes
test_genesymbol1 <- intersect(M4_string_db_ECC_merge_V717I_df$gene1symbol, M4_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M4/M4_string_db_genesymbol1_intersect_PPI_V717I.csv", row.names = F) 

#Note: checking which gene2symbols are M4 genes
test_genesymbol2 <- intersect(M4_string_db_ECC_merge_V717I_df$gene2symbol, M4_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M4/M4_string_db_genesymbol2_intersect_PPI_V717I.csv", row.names = F)

##############################################################################################
#4d. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 4 (APP vs. NDC)
##############################################################################################

#Included M4 genes in ECC + ReMap interactions + additional TF's with only
#M4 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_V717I <- cbind(gene2symbol = rownames(reslimma_V717I), reslimma_V717I)
rownames(reslimma_V717I) <- 1:nrow(reslimma_V717I)

#Subset M4 genes into only stringDB 
#Merge two data frames by ID
M4_string_db_ECC_merge_V717I_df <- merge(new_ECC_ReMap_M4,reslimma_V717I,by="gene2symbol")

#Note: checking which gene1symbols are M4 genes
test_genesymbol1 <- intersect(M4_string_db_ECC_merge_V717I_df$gene1symbol, M4_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M4/M4_string_db_genesymbol1_intersect_ECC_ReMap_V717I.csv", row.names = F) 

#Note: checking which gene2symbols are M4 genes
test_genesymbol2 <- intersect(M4_string_db_ECC_merge_V717I_df$gene2symbol, M4_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M4/M4_string_db_genesymbol2_intersect_ECC_ReMap_V717I.csv", row.names = F) 

#New Table with EDGES information - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M4_string_db_ECC_merge_V717I_df, file="~/M4/M4_ECC_ReMap_merge_V717I.csv", row.names = F) 

##################################################################
#4e. Enrichment for PPI interactions for Module 4 (PSEN1 vs. NDC)
##################################################################

#Included M4 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_A79V <- cbind(gene2symbol = rownames(reslimma_A79V), reslimma_A79V)
rownames(reslimma_A79V) <- 1:nrow(reslimma_A79V)

#Subset M4 genes into only stringDB 
#Merge two data frames by ID
M4_string_db_ECC_merge_A79V_df <- merge(new_stringdb_M4,reslimma_A79V,by="gene2symbol")


#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M4_string_db_ECC_merge_A79V_df, file="~/M4/M4_string_db_merge_A79V_v3.csv", row.names = F) 

#Note: checking which gene1symbols are M4 genes
test_genesymbol1_A79V <- intersect(M4_string_db_ECC_merge_A79V_df$gene1symbol, M4_genes$genesymbol1)  
View(test_genesymbol1_A79V)

#Convert character to data table 
test_genesymbol1_A79V <- as.data.table(test_genesymbol1_A79V)
write.csv(test_genesymbol1_A79V, file="~/M4/M4_string_db_genesymbol1_intersect_PPI_A79V.csv", row.names = F) 

#Note: checking which gene2symbols are M4 genes
test_genesymbol2_A79V <- intersect(M4_string_db_ECC_merge_A79V_df$gene2symbol, M4_genes$genesymbol1)  
View(test_genesymbol2_A79V)

#Convert character to data table 
test_genesymbol2_A79V <- as.data.table(test_genesymbol2_A79V)
write.csv(test_genesymbol2_A79V, file="~/M4/M4_string_db_genesymbol2_intersect_PPI_A79V.csv", row.names = F) 

###############################################################################################
#4f. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 4 (PSEN1 vs. NDC)
###############################################################################################

#Included M4 genes in ECC + ReMap interactions + additional TF's with only
#M4 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_A79V <- cbind(gene2symbol = rownames(reslimma_A79V), reslimma_A79V)
rownames(reslimma_A79V) <- 1:nrow(reslimma_A79V)

#Subset M4 genes into only stringDB 
#Merge two data frames by ID
M4_string_db_ECC_merge_A79V_df <- merge(new_ECC_ReMap_M4,reslimma_A79V,by="gene2symbol")


#Note: checking which gene1symbols are M4 genes
test_genesymbol1 <- intersect(M4_string_db_ECC_merge_A79V_df$gene1symbol, M4_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M4/M4_string_db_genesymbol1_intersect_ECC_ReMap_A79V.csv", row.names = F) 

#Note: checking which gene2symbols are M4 genes
test_genesymbol2 <- intersect(M4_string_db_ECC_merge_A79V_df$gene2symbol, M4_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M4/M4_string_db_genesymbol2_intersect_ECC_ReMap_A79V.csv", row.names = F) 

#New Table with EDGES information - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M4_string_db_ECC_merge_A79V_df, file="~/M4/M4_ECC_ReMap_merge_A79V.csv", row.names = F) 

######################################################################################

##################################################################
#4g. Enrichment for PPI interactions for Module 4 (PSEN2 vs. NDC)
##################################################################

#Included M4 genes in ONLY PPI interactions

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_N141I <- cbind(gene2symbol = rownames(reslimma_N141I), reslimma_N141I)
rownames(reslimma_N141I) <- 1:nrow(reslimma_N141I)

#Subset M4 genes into only stringDB 
#Merge two data frames by ID
M4_string_db_ECC_merge_N141I_df <- merge(new_stringdb_M4,reslimma_N141I,by="gene2symbol")

#New Table with EDGES information - stringDB PPI interactions
#Export as a .csv file
write.csv(M4_string_db_ECC_merge_N141I_df, file="~/M4/M4_string_db_merge_N141I_v3.csv", row.names = F) 

#Note: checking which gene1symbols are M4 genes
test_genesymbol1_N141I <- intersect(M4_string_db_ECC_merge_N141I_df$gene1symbol, M4_genes$genesymbol1)  
View(test_genesymbol1_N141I)

#Convert character to data table 
test_genesymbol1_N141I <- as.data.table(test_genesymbol1_N141I)
write.csv(test_genesymbol1_N141I, file="~/M4/M4_string_db_genesymbol1_intersect_PPI_N141I.csv", row.names = F) 

#Note: checking which gene2symbols are M4 genes
test_genesymbol2_N141I <- intersect(M4_string_db_ECC_merge_N141I_df$gene2symbol, M4_genes$genesymbol1)  
View(test_genesymbol2_N141I)

#Convert character to data table 
test_genesymbol2_N141I <- as.data.table(test_genesymbol2_N141I)
write.csv(test_genesymbol2_N141I, file="~/M4/M4_string_db_genesymbol2_intersect_PPI_N141I.csv", row.names = F) 


##############################################################################################
#4h. Enrichment for ENCODE/ChEA + ReMap and other TF interactions for Module 4 (PSEN2 vs. NDC)
##############################################################################################

#Included M4 genes in ECC + ReMap interactions + additional TF's with only
#M4 genes as target genes

#You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Move index column to first row
#Source: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
reslimma_N141I <- cbind(gene2symbol = rownames(reslimma_N141I), reslimma_N141I)
rownames(reslimma_N141I) <- 1:nrow(reslimma_N141I)

#Subset M4 genes into only stringDB 
#Merge two data frames by ID
M4_string_db_ECC_merge_N141I_df <- merge(new_ECC_ReMap_M4,reslimma_N141I,by="gene2symbol")


#Note: checking which gene1symbols are M4 genes
test_genesymbol1 <- intersect(M4_string_db_ECC_merge_N141I_df$gene1symbol, M4_genes$genesymbol1)  
View(test_genesymbol1)

#Convert character to data table 
test_genesymbol1 <- as.data.table(test_genesymbol1)
write.csv(test_genesymbol1, file="~/M4/M4_string_db_genesymbol1_intersect_ECC_ReMap_N141I.csv", row.names = F) 

#Note: checking which gene2symbols are M4 genes
test_genesymbol2 <- intersect(M4_string_db_ECC_merge_N141I_df$gene2symbol, M4_genes$genesymbol1)  
View(test_genesymbol2)

#Convert character to data table 
test_genesymbol2 <- as.data.table(test_genesymbol2)
write.csv(test_genesymbol2, file="~/M4/M4_string_db_genesymbol2_intersect_ECC_ReMap_N141I.csv", row.names = F) 

#New Table with EDGES information  - ECC + ReMap TF-gene interactions
#Export as a .csv file
write.csv(M4_string_db_ECC_merge_N141I_df, file="~/M4/M4_ECC_ReMap_merge_N141I.csv", row.names = F) 

#####################################################################
#5. Perform fGSEA enrichment for the different modules (M1, M3, + M4)
#####################################################################

M1_genes_nodes <- M1_genes
M3_genes_nodes <- M3_genes
M4_genes_nodes <- M4_genes

#Get the unique gene symbols from Module 1 
M1_genes_nodes_unique <- unique(M1_genes_nodes)

#Get the unique gene symbols from Module 3 
M3_genes_nodes_unique <- unique(M3_genes_nodes)

#Get the unique gene symbols from Module 4 
M4_genes_nodes_unique <- unique(M4_genes_nodes)

#Create new cem.modules object 
cem.modules <- rbind(M1_genes_nodes_unique, M3_genes_nodes_unique)
cem.modules <- rbind(cem.modules, M4_genes_nodes_unique)

#CEMiTool on filtered genes (includes both source and first neighbor genes)
modules <- unique(cem.modules[, 'modules'])

#Create the function for gene sets
gene_sets <- lapply(modules, function(mod){
  return(cem.modules[cem.modules[, 'modules']==mod, 'genes'])
})

names(gene_sets) <- modules
cem.modules.gmt <- gene_sets

##########################################################################
#5a. Perform topTable for V717I vs. NDC to get ranked gene list
##########################################################################

reslimma_V717I <- topTable(efit, adjust.method="BH", sort.by = "t", coef = 1, n = Inf)

#Convert your limma topTable results to a data table format
reslimma_V717I_dt <- reslimma_V717I

#Move the first column as a column called SYMBOL
reslimma_V717I_dt <- cbind(SYMBOL = rownames(reslimma_V717I_dt), reslimma_V717I_dt)
rownames(reslimma_V717I_dt) <- 1:nrow(reslimma_V717I_dt)

#Remove rows with an empty or duplicated Symbol
reslimma_V717I_dt <- subset(reslimma_V717I_dt, SYMBOL != "" )
reslimma_V717I_dt <- subset(reslimma_V717I_dt, ! duplicated(SYMBOL))

#Create data table
reslimma_V717I_dt <- data.table(reslimma_V717I_dt)

#Tidy up and order by t metric
ranks_limma_V717I <- as_tibble(reslimma_V717I_dt[order(t), list(SYMBOL, t)])

#Deframe the data table
ranks_t_V717I <- deframe(ranks_limma_V717I)

##########################################################################
#5b. Perform topTable for A79V vs. NDC to get ranked gene list
##########################################################################

reslimma_A79V <- topTable(efit, adjust.method="BH", sort.by = "t", coef = 2, n = Inf)

#Convert your limma topTable results to a data table format
reslimma_A79V_dt <- reslimma_A79V

#Move the first column as a column called SYMBOL
reslimma_A79V_dt <- cbind(SYMBOL = rownames(reslimma_A79V_dt), reslimma_A79V_dt)
rownames(reslimma_A79V_dt) <- 1:nrow(reslimma_A79V_dt)

#Remove rows with an empty or duplicated Symbol
reslimma_A79V_dt <- subset(reslimma_A79V_dt, SYMBOL != "" )
reslimma_A79V_dt <- subset(reslimma_A79V_dt, ! duplicated(SYMBOL))

#Create data table
reslimma_A79V_dt <- data.table(reslimma_A79V_dt)

#Tidy up and order by t metric
ranks_limma_A79V <- as_tibble(reslimma_A79V_dt[order(t), list(SYMBOL, t)])

#Deframe the data table
ranks_t_A79V <- deframe(ranks_limma_A79V)

###########################################################################
#5c. Perform topTable for N141I vs. NDC to get ranked gene list
###########################################################################

reslimma_N141I <- topTable(efit, adjust.method="BH", sort.by = "t", coef = 3, n = Inf)

#Convert your limma topTable results to a data table format
reslimma_N141I_dt <- reslimma_N141I

#Move the first column as a column called SYMBOL
reslimma_N141I_dt <- cbind(SYMBOL = rownames(reslimma_N141I_dt), reslimma_N141I_dt)
rownames(reslimma_N141I_dt) <- 1:nrow(reslimma_N141I_dt)

#Remove rows with an empty or duplicated Symbol
reslimma_N141I_dt <- subset(reslimma_N141I_dt, SYMBOL != "" )
reslimma_N141I_dt <- subset(reslimma_N141I_dt, ! duplicated(SYMBOL))

#Create data table
reslimma_N141I_dt <- data.table(reslimma_N141I_dt)

#Tidy up and order by t metric
ranks_limma_N141I <- as_tibble(reslimma_N141I_dt[order(t), list(SYMBOL, t)])

#Deframe the data table
ranks_t_N141I <- deframe(ranks_limma_N141I)

#####################################
#5d. Loading the Gene Set Databases
#####################################

#Get modules for only module genes in the modules
cem.modules_fgsea <- module_genes(cem)
modules <- unique(cem@module[, 'modules'])
gene_sets_fgsea <- lapply(modules, function(mod){
  return(cem.modules_fgsea[cem.modules_fgsea[, 'modules']==mod, 'genes'])
})
#
cem.modules_fgsea <- na.omit(cem.modules_fgsea)
names(gene_sets_fgsea) <- modules
cem.modules_fgsea.gmt <- gene_sets_fgsea

################################################################
#5e. fGSEA with t-value Limma Ranking for APP-V717I vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_V717I as a ranking for fGSEA

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes_V717I.cem <- fgsea::fgseaMultilevel(pathways=cem.modules_fgsea.gmt, stats=ranks_t_V717I, absEps = 0)

# Make the Results tidy
fgseaRestidy_V717I.cem <- fgseaRes_V717I.cem %>%
  as_tibble() %>%
  arrange(desc(NES))

#Save results to .csv file
write.csv(fgseaRestidy_V717I.cem, file = "~/Results_fGSEA_modules_V717I_NDC.csv", row.names = TRUE)

################################################################
#5f. fGSEA with t-value Limma Ranking for PSEN1-A79V vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_A79V as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes_A79V.cem <- fgsea::fgseaMultilevel(pathways=cem.modules_fgsea.gmt, stats=ranks_t_A79V, absEps = 0)

# Make the Results tidy
fgseaRes_A79V.cem <- fgseaRes_A79V.cem %>%
  as_tibble() %>%
  arrange(desc(NES))

#Save results to .csv file
write.csv(fgseaRes_A79V.cem, file = "~/Results_fGSEA_modules_A79V_NDC.csv", row.names = TRUE)

################################################################
#5g. fGSEA with t-value Limma Ranking for PSEN2-N141I vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_N141I as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes_N141I.cem <- fgsea::fgseaMultilevel(pathways=cem.modules_fgsea.gmt, stats=ranks_t_N141I, absEps = 0)

# Make the Results tidy
fgseaRes_N141I.cem <- fgseaRes_N141I.cem %>%
  as_tibble() %>%
  arrange(desc(NES))

#Save results to .csv file
write.csv(fgseaRes_N141I.cem, file = "~/Results_fGSEA_modules_N141I_NDC.csv", row.names = TRUE)

#########################################################
#6. Run hypergeometric test for modules (M1, M3 and M4)
#########################################################

#Import the latest MSigDB XML File
msig <- tmodImportMSigDB(file.path("~/msigdb", "msigdb_v7.4.xml"))

Hallmark.sel <- msig$MODULES$Category == "H"
GOBP.sel <- msig$MODULES$Subcategory == "C5_GO:BP"
Reactome.sel <- msig$MODULES$Subcategory == "CP:REACTOME"
GOBP.exp.sel <- msig$MODULES$Subcategory == "GO:BP"
msetECC
msetReMap

##########################################
#6a. Run hypergeometric test for Module 1
##########################################

mod_genes <- cem.modules.gmt$M1

tmodHG_GOBP.M1 <- tmodHGtest(mod_genes, all_genes, mset = msig[GOBP.sel])
head(tmodHG_GOBP.M1,50)
write.csv(tmodHG_GOBP.M1, file = "~/M1/tmodHG_GOBP.M1.csv") 

tmodHG_Hallmark.M1 <- tmodHGtest(mod_genes, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark.M1,50)
write.csv(tmodHG_Hallmark.M1, file = "~/M1/tmodHG_Hallmark.M1.csv") 

tmodHG_ECC.M1 <- tmodHGtest(mod_genes, all_genes, mset = msetECC)
head(tmodHG_ECC.M1,50)
write.csv(tmodHG_ECC.M1, file = "~/M1/tmodHG_ECC.M1.csv") 

tmodHG_ReMap.M1 <- tmodHGtest(mod_genes, all_genes, mset = msetReMap)
head(tmodHG_ReMap.M1,50)
write.csv(tmodHG_ReMap.M1, file = "~/M1/tmodHG_ReMap.M1.csv") 

##########################################
#6b. Run hypergeometric test for Module 3
##########################################

#mod_genes <- cem.modules.gmt$M3
mod_genes_M3 <- cem.modules.gmt$M3

tmodHG_GOBP.M3 <- tmodHGtest(mod_genes_M3, all_genes, mset = msig[GOBP.sel])
head(tmodHG_GOBP.M3,50)
write.csv(tmodHG_GOBP.M3, file = "~/M3/tmodHG_GOBP.M3.csv") 

tmodHG_Hallmark.M3 <- tmodHGtest(mod_genes_M3, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark.M3,50)
write.csv(tmodHG_Hallmark.M3, file = "~/M3/tmodHG_Hallmark.M3.csv") 

tmodHG_ECC.M3 <- tmodHGtest(mod_genes_M3, all_genes, mset = msetECC)
head(tmodHG_ECC.M3,50)
write.csv(tmodHG_ECC.M3, file = "~/M3/tmodHG_ECC.M3.csv") 

tmodHG_ReMap.M3 <- tmodHGtest(mod_genes_M3, all_genes, mset = msetReMap)
head(tmodHG_ReMap.M3,50)
write.csv(tmodHG_ReMap.M3, file = "~/M3/tmodHG_ReMap.M3.csv") 

##########################################
#6c. Run hypergeometric test for Module 4
##########################################

#mod_genes <- cem.modules.gmt$M3
mod_genes_M4 <- cem.modules.gmt$M4

tmodHG_GOBP.M4 <- tmodHGtest(mod_genes_M4, all_genes, mset = msig[GOBP.sel])
head(tmodHG_GOBP.M4,50)
write.csv(tmodHG_GOBP.M4, file = "~/M4/tmodHG_GOBP.M4.csv") 

tmodHG_Reactome.M4 <- tmodHGtest(mod_genes_M4, all_genes, mset = msig[Reactome.sel])
head(tmodHG_Reactome.M4,50)
write.csv(tmodHG_Reactome.M4, file = "~/M4/tmodHG_Reactome.M4.csv") 

tmodHG_Hallmark.M4 <- tmodHGtest(mod_genes_M4, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark.M4,50)
write.csv(tmodHG_Hallmark.M4, file = "~/M4/tmodHG_Hallmark.M4.csv") 

tmodHG_KEGG.M4 <- tmodHGtest(mod_genes_M4, all_genes, mset = msig[KEGG.sel])
head(tmodHG_KEGG.M4,50)
write.csv(tmodHG_KEGG.M4, file = "~/M4/tmodHG_KEGG.M4.csv") 

tmodHG_ECC.M4 <- tmodHGtest(mod_genes_M4, all_genes, mset = msetECC)
head(tmodHG_ECC.M4,50)
write.csv(tmodHG_ECC.M4, file = "~/M4/tmodHG_ECC.M4.csv") 

tmodHG_ReMap.M4 <- tmodHGtest(mod_genes_M4, all_genes, mset = msetReMap)
head(tmodHG_ReMap.M4,50)
write.csv(tmodHG_ReMap.M4, file = "~/M4/tmodHG_ReMap.M4.csv") 

##########################################
