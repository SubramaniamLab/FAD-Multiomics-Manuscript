#RNA-Seq PCA and tSNE analysis 
#Valdez et. al 2023 Molecular Psychiatry Submission
#Load packages
library("edgeR")
library("limma")
library("rhdf5")
library("readr")
library("rjson")
library("tximport")
library("ggplot2")
library("ensembldb")
library("org.Hs.eg.db")
library("tibble")
library("qusage")
library("tidyverse")
library("fgsea")
library("RColorBrewer")
library("Glimma")
library("VennDiagram")
library("statmod")
library("tidyr")
library("dplyr")
library("httr")
library("jsonlite")
library("biomaRt")
library("stats")
library("ggplot2")
library("ggrepel")
library("tmod")
library("data.table")
library("Rtsne")

load("~/PSEN1.PSEN2.APP.txi_lsTPM.Rdata")
load("~/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")
######################################
#Get Raw, Unfiltered RNA-seq Count Data
########################################

#Get raw counts from the tximport object
txi_counts <- txi_lsTPM$counts

#Export raw count data ran by PV from R into .csv files
write.csv(txi_counts, '~/PSEN1_PSEN2_APP_NDC_samples_raw_unfiltered_counts.csv', append =FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

####################################################
#Creating PCA plot for Filtered data
####################################################

#Primer on PCA here:
#http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html
#https://rpubs.com/Mentors_Ubiqum/Transpose_Dataframe

#Transpose the counts matrix
#transpose the filtered data to have variables (genes) as columns and (samples) as rows
data_for_PCA_filtered <- as.data.frame((as.matrix(t(filteredExpr))))

#Dimensions should be [22310 genes] x [12 samples] 
dim(data_for_PCA_filtered)

#Can also do the following approach using R
#Source: https://cmdlinetips.com/2019/04/introduction-to-pca-with-r-using-prcomp/

#Transpose the PCA data matrix first
data_for_PCA_filt <- as.data.frame(t(as.matrix(data_for_PCA_filtered)))

#Compute PCA for filtered data
myPrcomp_filt <- prcomp(data_for_PCA_filt, scale. = TRUE)
summary(myPrcomp_filt)
autoplot(myPrcomp_filt)

#Compute the variance explained
var_explained_filt <- myPrcomp_filt$sdev^2/sum(myPrcomp_filt$sdev^2)
var_explained_filt[1:2]

#Make PCA scatter plot using PC1 and PC2 grouped by condition
#and replicate ID for filtered data

#Create column vector showing replicates
replicate_ID <- c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3")

#Create column vector showing conditions
condition <- c(rep("NDC",3), rep("V717I",3), rep("A79V",3), rep("N141I",3))

#Make PCA scatter plot using PC1 and PC2 grouped by diagnostic code
#and replicate ID for unfiltered data

#Increase the max overlaps for ggrepel
options(ggrepel.max.overlaps = 15)

#myPrcomp_filt_t <- t(myPrcomp_filt)

#Make PCA scatter plot using PC1 and PC2 grouped by condition
#and replicate ID for unfiltered data
myPrcomp_filt$x %>% 
  as.data.frame %>%
  rownames_to_column("SampleID") %>%
  add_column(subcondition = replicate_ID) %>%
  add_column(condition = condition) %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=condition),size=4) +
  geom_text_repel(aes(label =subcondition),
                  force = 25,
                  segment.size  = 0.2,
                  segment.color = "grey50") +
  scale_color_manual(values = c("#DC267F", "#785EF0", "#999999", "#FE6100")) +  
  labs(x=paste0("PC1: ",round(var_explained_filt[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained_filt[2]*100,1),"%")) +
  ggtitle ("PCA of Samples Based on Filtered Data", ) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24, face ="bold"), axis.text = element_text(size=20), axis.title = element_text(size=20), legend.position="bottom") #center title and bold face

################################################
#Retreive PC1 and PC2 Values to input for Prism
################################################

#Get values
#Get PCA scores
PCA_scores_filtered <- myPrcomp_filt$x

#Write the .csv file with all PCA scores 
write.csv(PCA_scores_filtered, '~/RNA-seq_PCA_scores_PSEN1_PSEN2_APP_NDC.csv', append =FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

################################################
#Create Venn overlap figures with nVennR
################################################

#Do topTable() for the APP-V717I vs. NDC 
V717I_DEG <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf) 
V717I_DEG <- V717I_DEG %>%
  dplyr::filter(adj.P.Val < 0.05)
V717I <- list(V717I = as.character(unique(V717I_DEG$SYMBOL)))

#Do topTable() for the PSEN1-A79V vs. NDC 
A79V_DEG <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf) 
A79V_DEG <- A79V_DEG %>%
  dplyr::filter(adj.P.Val < 0.05)
A79V <- list(A79V = as.character(unique(A79V_DEG$SYMBOL)))

#Do topTable() for the PSEN2-N141I vs. NDC  
N141I_DEG <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf) 
N141I_DEG <- N141I_DEG %>%
  dplyr::filter(adj.P.Val < 0.05)
N141I <- list(N141I = as.character(unique(N141I_DEG$SYMBOL)))

#Plot the nVenn Diagram
myNV3 <- plotVenn(list(V717I, A79V, N141I), sNames=c("V717I", "A79V", "N141I"),showPlot = T,nCycles = 1000000)
showSVG(myNV3, opacity=0.3,outFile = "~/All3_PSEN1_PSEN2_APP_DEG_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))

################################################
#Create raw counts table for input into diffTF 
################################################

#Create filtered raw counts data table
final_counts = as.data.frame(y$counts)

#Move the index columns into the first column called ENSEMBL
final_counts <- cbind(ENSEMBL = rownames(final_counts), final_counts)
View(final_counts)

#Source: https://stackoverflow.com/questions/6081439/changing-column-names-of-a-data-frame

#Rename names of colummns
colnames(final_counts) <- c("ENSEMBL","NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3",
                            "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3",
                            "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")

final_counts = as.data.frame(final_counts)

#Export raw unfiltered count data into .csv files
write.csv(final_counts, '~/PSEN1_PSEN2_APP_unfiltered_unnormalized_countdata.csv', append =FALSE, sep ="\t", dec = ".", row.names= FALSE, col.names = TRUE)

#######################################
#tSNE for Filtered log2 Normalized Data
#######################################

####################
#Filtering Genes
####################

# Filtering lowly expressed genes
# This function is used in tximport manual to remove lowly expressed genes
# Only kept genes with 5 counts/gene
keep = filterByExpr(y, min.count = 5)
y_filtered <- y[keep,]

#Obtain the column of gene symbols
final_genes = as.data.frame(y_filtered$genes["SYMBOL"])
final_genes2 = as.data.frame(y_filtered$genes)

#Convert to DGEList
y_filtered <- DGEList(y_filtered)

#Number of genes left after filtering (28,305 genes x 18 samples)
dim(y_filtered)

#Create filtered raw counts data table
final_counts_filtered = as.data.frame(y_filtered$counts)

#Move the index columns into the first column called ENSEMBL
final_counts_filtered <- cbind(ENSEMBL = rownames(final_counts_filtered), final_counts_filtered)
View(final_counts_filtered)

#Source: https://stackoverflow.com/questions/6081439/changing-column-names-of-a-data-frame

#Rename names of colummns
colnames(final_counts_filtered) <- c("ENSEMBL","NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3",
                                     "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3",
                                     "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")

final_counts_filtered = as.data.frame(final_counts_filtered)

#Create a new filtered object with gene symbols for tSNE

#Obtain filtered count matrix with gene symbols
#Combine the gene names and the integer expresson values together
final_counts_filtered_genes <- cbind(final_genes,filteredExpr2)
expression_data <- final_counts_filtered_genes

#Conduct tSNE analysis for filtered, log2 normalized data 
# x is normalized data (rows are genes, columns are samples)
# you may have to play with some of the parameters, e.g., perplexity

#Run tSNE using filtered, normalized data
tsne_out <- Rtsne(t(expression_data), pca=TRUE, perplexity = 5, theta = 0.1, 
                  verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                  max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                  num_threads = 1); 

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")

#Create plot with tSNE output
plot(tsne_out$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, xlab = "TSNE1", ylab = "TSNE2")

################################################
#tSNE Data Formatting For Four, Main Endotypes
################################################

#Read the .csv file with the list of genes for each endotype
custom_endo_data = read.csv("~/CustomEndotypes.csv")
cell_cycle_endo = data.frame(custom_endo_data$Cell_Cycle_Custom_Endtotype)
dediff_endo = data.frame(custom_endo_data$Dedifferentiation_Custom_Endotype)
neuron_lineage_endo = data.frame(custom_endo_data$Neuron_Lineage_Custom_Endtype)
neuron_function_endo = data.frame(custom_endo_data$Neuron_Function_Custom_Endtype)
inflammation_endo = data.frame(custom_endo_data$Inflammation)
G0_G1_endo = data.frame(custom_endo_data$G0_G1)
G1_S_endo = data.frame(custom_endo_data$G1_S)
G2_M_endo = data.frame(custom_endo_data$G2_M)

#Rename column
#Rename names of colummns
colnames(cell_cycle_endo) <- c("SYMBOL")
colnames(dediff_endo) <- c("SYMBOL")
colnames(neuron_lineage_endo) <- c("SYMBOL")
colnames(neuron_function_endo) <- c("SYMBOL")
colnames(inflammation_endo) <- c("SYMBOL")
colnames(G0_G1_endo) <- c("SYMBOL")
colnames(G1_S_endo) <- c("SYMBOL")
colnames(G2_M_endo) <- c("SYMBOL")

###############################
#Subsetting Cell Cycle Counts
###############################

#Subset counts into only cell cycle genes
#Merge two data frames by ID
#883 genes x 18 samples
cell_cycle_merge <- merge(expression_data_genes,cell_cycle_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
cell_cycle_merge <- as.data.frame(cell_cycle_merge[,-1], row.names = cell_cycle_merge[,1])

write.csv(cell_cycle_merge, '~/data_for_PCA_filtered_Cell_Cycle_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

######################################
#Subsetting Dedifferentiation Counts
######################################

#Subset counts into only dedifferentiation genes
#Merge two data frames by ID
#482 genes x 18 samples
dediff_merge <- merge(expression_data_genes,dediff_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
dediff_merge <- as.data.frame(dediff_merge[,-1], row.names = dediff_merge[,1])

write.csv(dediff_merge, '~/data_for_PCA_filtered_Dedifferentiation_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

######################################
#Subsetting Neuron Lineage Counts
######################################

#Subset counts into only neuron lineage genes
#Merge two data frames by ID
#1165 genes x 18 samples
neuron_lineage_merge <- merge(expression_data_genes,neuron_lineage_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
neuron_lineage_merge <- as.data.frame(neuron_lineage_merge[,-1], row.names = neuron_lineage_merge[,1])

write.csv(neuron_lineage_merge, '~/data_for_PCA_filtered_Neuron_Lineage_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

######################################
#Subsetting Neuron Function Counts
######################################

#Subset counts into only neuron function genes
#Merge two data frames by ID
#513 genes x 18 samples
neuron_function_merge <- merge(expression_data_genes,neuron_function_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
neuron_function_merge <- as.data.frame(neuron_function_merge[,-1], row.names = neuron_function_merge[,1])

write.csv(neuron_function_merge, '~/data_for_PCA_filtered_Neuron_Function_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

######################################
#Subsetting Inflammation Counts
######################################

#Subset counts into only neuron function genes
#Merge two data frames by ID
#532 genes x 18 samples
inflammation_merge <- merge(expression_data_genes,inflammation_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
inflammation_merge <- as.data.frame(inflammation_merge[,-1], row.names = inflammation_merge[,1])

write.csv(inflammation_merge, '~/data_for_PCA_filtered_Inflammation_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

######################################
#Subsetting G0_G1 Counts
######################################

#Subset counts into only neuron function genes
#Merge two data frames by ID
#39 genes x 18 samples
G0_G1_merge <- merge(expression_data_genes,G0_G1_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
G0_G1_merge <- as.data.frame(G0_G1_merge[,-1], row.names = G0_G1_merge[,1])

write.csv(G0_G1_merge, ~'/data_for_PCA_filtered_G0_G1_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

######################################
#Subsetting G1_S Counts
######################################

#Subset counts into only neuron function genes
#Merge two data frames by ID
#258 genes x 18 samples
G1_S_merge <- merge(expression_data_genes,G1_S_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
G1_S_merge <- as.data.frame(G1_S_merge[,-1], row.names = G1_S_merge[,1])

write.csv(G1_S_merge, '~/data_for_PCA_filtered_G1_S_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

######################################
#Subsetting G2_M Counts
######################################

#Subset counts into only neuron function genes
#Merge two data frames by ID
#273 genes x 18 samples
G2_M_merge <- merge(expression_data_genes,G2_M_endo,by="SYMBOL")

#Move the first column as index row names
#Source: https://stackoverflow.com/questions/45526629/convert-first-column-in-data-frame-to-row-index
G2_M_merge <- as.data.frame(G2_M_merge[,-1], row.names = G2_M_merge[,1])

write.csv(G2_M_merge, '~/data_for_PCA_filtered_G2_M_PV.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

##################
#Cell Cycle tSNE
##################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_cell_cycle <- Rtsne(t(cell_cycle_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                             verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                             max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                             num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")

# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_cell_cycle$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "Cell Cycle Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()

#########################
#Dedifferentiation tSNE
#########################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_dediff <- Rtsne(t(dediff_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                         verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                         max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                         num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")


# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_dediff$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "Dedifferentiation Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()

#######################
#Neuron Lineage tSNE
#######################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_neuron_lineage <- Rtsne(t(neuron_lineage_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                                 verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                                 max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                                 num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")


# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_neuron_lineage$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "Neuron Lineage Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()

#######################
#Neuron Function tSNE
#######################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_neuron_function <- Rtsne(t(neuron_function_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                                  verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                                  max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                                  num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")


# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_neuron_function$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "Neuron Function Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()

save.image("~/Desktop/tSNE_2021/tSNE_tximport_filtered_v3.RData")

#######################
#Inflammation tSNE
#######################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_inflammation <- Rtsne(t(inflammation_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                               verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                               max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                               num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")


# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_inflammation$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "Inflammation Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()

#######################
#G0_G1 tSNE
#######################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_G0_G1 <- Rtsne(t(G0_G1_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                        verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                        max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                        num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")


# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_G0_G1$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "G0_G1 Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()


#######################
#G1_S tSNE
#######################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_G1_S <- Rtsne(t(G1_S_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                       verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                       max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                       num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")


# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_G1_S$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "G1_S Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()

#######################
#G2_M tSNE
#######################

#Run tSNE using filtered, log2 nnormalized data
tsne_out_G2_M <- Rtsne(t(G2_M_merge), pca=TRUE, perplexity = 5, theta = 0.0, 
                       verbose = TRUE, partial_pca = TRUE, initial_dims = 17,  
                       max_iter = 5000, pca_center = TRUE, pca_scale = TRUE, 
                       num_threads = 1);

#For sample colors, create a vector of colors (as many elements as number of samples) and use col argument in the plot function.
cl = c(rep("black",3),rep("orange",3),  rep("red", 3), rep("blue", 3), rep("purple", 3), rep("green", 3))

#sample_labels contains the labels i.e., class for the samples
sample_labels = c("NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3", "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3", "H163R_1", "H163R_2", "H163R_3", "A431E_1", "A431E_2", "A431E_3")


# Add extra space to right of plot area; change clipping to figure
par(mar=c(4, 4, 3, 10), xpd=FALSE)  

#Create plot with tSNE output
plot(tsne_out_G2_M$Y,col=cl[as.factor(sample_labels)], asp=1, pch = 16, main = "G2_M Endotype", xlab = "TSNE1", ylab = "TSNE2")

# Add a legend
legend("topright", inset=c(-0.35,0), legend = c("NDC", "V717I", "A79V", "N141I", "H163R", "A431E") , 
       col = c("purple", "green", "orange", "blue", "red", "black") , bty = "n", xpd=TRUE,mar(c(-5,1,1,1)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

dev.off()

#######################
#Save the .RData file
#######################

save.image('~/PSEN1_PSEN2_APP_RNAseq_Analysis_tSNE.RData')