#RNA-Seq limma-voom DGE Analysis pipeline 
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
library("data.table")
library("RColorBrewer")
library("Glimma")
library("VennDiagram")
library("statmod")
library("tidyr")
library("httr")
library("jsonlite")
library("biomaRt")
library("stats")
library("ggrepel")

############################################################
#Ensembl Transcript to Gene ID's Conversion
############################################################
#--------------------------------------------------------------------------------------------------------
# Get the transcript to Gene IDs for Ensembl
# Use bioMart - make sure the host matches the version you use, which here is GRCh38.99
#useMart enables connection to a specified BioMart database
#Ensembl version 99 
martGRCh38.99 = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                 dataset = "hsapiens_gene_ensembl",
                                 host = 'jan2020.archive.ensembl.org',
                                 path="/biomart/martservice")
#getBM is the main biomaRt query function and this usually retrieves the user's specified attributes
#from the BioMart database
GRCh38.99t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id_version",
                                             "ensembl_gene_id"),
                              mart = martGRCh38.99)
#dplyr::rename function is used to rename columns
GRCh38.99t2g = dplyr::rename(GRCh38.99t2g, 
                             TXNAME = ensembl_transcript_id_version, 
                             ENSEMBL = ensembl_gene_id)
head(GRCh38.99t2g)
# save Ensembl mart object
save(martGRCh38.99, "~/gmtfiles/martGRCh38.99.Rdata")
############################################################
#Accessing Files Needed in Base Directory
############################################################
#--------------------------------------------------------------------------------------------------------
# Set the base directory containing kallisto output files
base_dir <- "~/kallisto_dir"
#Get samples for  model from 3 column txt file (sample name, condition, replicate)
samples <- read.table(file.path(base_dir, "PSEN1.PSEN2.APP_info.txt"), header = TRUE, stringsAsFactors=FALSE)
#For Kallisto, describe the path to find the quant.sf files
files <- file.path(base_dir, samples$run_sample, "abundance.h5")
#Apply the sample names to "files"
names(files) <- paste0(c(samples$run_sample))
# Check if all files exist
all(file.exists(files))
############################################################
#Import Counts Using Tximport
############################################################
#--------------------------------------------------------------------------------------------------------
#Import the abundance/counts measurements using tximport

txi_lsTPM = tximport(files, 
                     type = "kallisto", 
                     tx2gene = GRCh38.99t2g, 
                     countsFromAbundance = "lengthScaledTPM")

# Check the head of txi_lsTPM (TPM = transcripts per million)
head(txi_lsTPM$counts)
names(txi_lsTPM)

# Create a sample table with the list of conditions
sampleTable <- data.frame(condition=factor(rep(c(samples$condition))))

dim(txi_lsTPM$counts)

save(txi_lsTPM, file = "~/PSEN1.PSEN2.APP.txi_lsTPM.Rdata")
############################################################
#Limma-Voom Differential Expression
############################################################

#--------------------------------------------------------------------------------------------------------
## PERFORM DIFFERENTIAL EXPRESSION WITH LIMMA-VOOM ##
# Convert counts to DGEList
y <- DGEList(txi_lsTPM$counts,
             lib.size = colSums(txi_lsTPM$counts),
             norm.factors = calcNormFactors(txi_lsTPM$counts),
             samples = samples$sample,
             group = samples$condition)

#--------------------------------------------------------------------------------------------------------
#Unloading dplyr package before using the select() method
detach("package:dplyr", unload=TRUE)

#Create a Homo Sapiens annotation from the org.Hs.eg.db database 
Hs_ann = select(org.Hs.eg.db,
                keys=rownames(y$counts),
                columns=c("ENTREZID","SYMBOL"),
                keytype="ENSEMBL",
                multiVals="first")

# Remove duplicated terms
Hs_ann <- Hs_ann[!duplicated(Hs_ann[,1]),]
head(Hs_ann)
dim(Hs_ann)

#Apply the annotation to your limma object "y"
#Match gene symbols to Ensembl ID's
y$genes <- Hs_ann

# View the library size for each sample
y$samples

#Number of genes (66023 12)
dim(y)

#--------------------------------------------------------------------------------------------------------
# Load a color palette of 50 colors to be used for plots
myPalette <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))

#Convert counts to cpm and log
unfilteredExpr <- cpm(y, log=T)

#Transpose the unfiltered expression matrix
unfilteredExpr <- as.data.frame(t(as.matrix(unfilteredExpr)))
View(unfilteredExpr)

#Rename the columns
(setattr(unfilteredExpr, "row.names", c("NDC_1", "NDC_2", "NDC_3",
                                        "V717I_1", "V717I_2", "V717I_3",
                                        "A79V_1", "A79V_2", "A79V_3",
                                        "N141I_1", "N141I_2", "N141I_3")))


#Create a new unfiltered Expression object that transposes the PCA data matrix
unfilteredExpr2 <- t(unfilteredExpr)

#Export unfiltered count data ran by PV from R into .txt files
write.csv(unfilteredExpr2, 'unfiltered_countdata.csv', append =FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)


# Plot the density of unfiltered gene expression for all samples within groups
plotDensities(unfilteredExpr2, group=samples$condition, col=myPalette[1:4], legend="topright", main ="Distribution by Conditions of Unfiltered Data")

#--------------------------------------------------------------------------------------------------------
# Add extra space to right of plot area; change clipping to figure
#(bottom, left, top, right)
par(mar=c(4, 2, 3, 8), xpd=FALSE)  

# Prepare a vector of colors with specific color for NDC samples, V717I_alt samples, A79V_alt samples and A141I samples
colors = c(rep(rgb(0.1,0.1,0.7,0.5),3),rep("#6BB2FF",3),  rep("#00CED1", 3), rep("#BC8F8F", 3))

#Change size of the labels for the x-axis
par(cex.axis=0.6) # is for x-axis

#Make boxplots to compare for unnormalized data?plot
boxplot(unfilteredExpr2, names = colnames(unfilteredExpr2), col = colors, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalized logCPM")

## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(unfilteredExpr2),col="blue")

# Add a legend
legend("topright", inset=c(-0.3,0), legend = c("NDC", "V717I", "A79V", "N141I") , 
       col = c(rgb(0.1,0.1,0.7,0.5) , "#6BB2FF", "#00CED1", "#BC8F8F") , bty = "n", xpd=TRUE,mar(c(65,20,10,25)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

#Return to default mode
dev.off()

#--------------------------------------------------------------------------------------------------------
# Filtering lowly expressed genes with edgeR
keep = filterByExpr(y)
y <- y[keep,]
y <- DGEList(y)

#Number of genes left after filtering (22310 12)
dim(y)

#--------------------------------------------------------------------------------------------------------
# Calculating normalization factors
y <- calcNormFactors(y)
#Add genes to the 'y' object
y$genes <- Hs_ann
#Add samples to the 'y' object
y$samples

save(y, file = "PSEN1.PSEN2.APP.y.Rdata")

#---------------------------------------------------------------------------------------------------
#Plot the density of filtered gene expression for all samples within groups
filteredExpr <- cpm(y, log=T)

#Transpose the filtered expression matrix
filteredExpr <- as.data.frame(t(as.matrix(filteredExpr)))
View(filteredExpr)

#Rename the columns
(setattr(filteredExpr, "row.names", c("NDC_1", "NDC_2", "NDC_3",
                                      "V717I_1", "V717I_2", "V717I_3",
                                      "A79V_1", "A79V_2", "A79V_3",
                                      "N141I_1", "N141I_2", "N141I_3")))

#Create a new unfiltered Expression object that transposes the PCA data matrix
filteredExpr2 <- t(filteredExpr)

#Export filtered count data ran by PV from R into .txt files
write.csv(filteredExpr2, 'filtered_countdata.csv', append =FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

#--------------------------------------------------------------------------------------------------------
#Add extra space to right of plot area; change clipping to figure
#(bottom, left, top, right)
par(mar=c(4, 2, 3, 8), xpd=FALSE) 

#Plot the density of unfiltered gene expression for all samples within groups
plotDensities(filteredExpr2, group=samples$condition, col=myPalette[1:4],, legend = "topright", main ="Distribution by Conditions of Filtered Data")

#Prepare a vector of colors with specific color for NDC samples, V717I, V717I_alt samples, A79V, A79V_alt samples and A141I, A141I_alt samples
colors = c(rep(rgb(0.1,0.1,0.7,0.5),3),rep("#6BB2FF",3),  rep("#00CED1", 3), rep("#BC8F8F", 3))

#Change size of the labels for the x-axis
par(cex.axis=0.6) # is for x-axis

#Make boxplots to compare for unnormalized data?plot
boxplot(filteredExpr2, names = colnames(filteredExpr2), col = colors, xlab="", ylab="Log2 counts per million",las=2,main="Normalized logCPM")

## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(filteredExpr2),col="blue")

# Add a legend
legend("topright", inset=c(-0.3,0), legend = c("NDC", "V717I", "A79V", "N141I") , 
       col = c(rgb(0.1,0.1,0.7,0.5) , "#6BB2FF", "#00CED1", "#BC8F8F") , bty = "n", xpd=TRUE,mar(c(65,20,10,25)), pch=20 , pt.cex = 3, cex = 1, horiz = FALSE)

#Return to default mode
dev.off()

#--------------------------------------------------------------------------------------------------------
#create the group and design 
#Condition1 is the NDC or the non-standard control
#Condition2 represents the V717I mutation from APP
#Condition3 represents the A79V mutation from PSEN1
#Condition4 represents the N141I mutation from PSEN2
group <- factor(samples$condition,levels=c("NDC","V717I", "A79V", "N141I"))

#Multidimensional scaling (MDS) plot
plotMDS(y, col = as.numeric(group))

#--------------------------------------------------------------------------------------------------------
##How to make a PCA plot for unfiltered data##

#Resources for PCA here:
#http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html
#https://rpubs.com/Mentors_Ubiqum/Transpose_Dataframe

#Transpose the counts matrix
# transpose the filtered data to have variables (genes) as columns and (samples) as rows
data_for_PCA_unfiltered <- as.data.frame((as.matrix(unfilteredExpr2)))

#Dimensions should be [66023 genes] x [12 samples]
dim(data_for_PCA_unfiltered)

#Write the data matrix into a text file
write.csv(data_for_PCA_unfiltered, 'data_for_PCA_unfiltered.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

#--------------------------------------------------------------------------------------------------------
##How to make a PCA plot for filtered data##

#Resources for PCA here:
#http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html
#https://rpubs.com/Mentors_Ubiqum/Transpose_Dataframe

#Transpose the counts matrix
# transpose the filtered data to have variables (genes) as columns and (samples) as rows
data_for_PCA_filtered <- as.data.frame((as.matrix(filteredExpr2)))

#Dimensions should be [22310 genes] x [12 samples]
dim(data_for_PCA_filtered)

#Write the data matrix into a text file
write.csv(data_for_PCA_filtered, 'data_for_PCA_filtered.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

#Conduct PCA analysis using Python to make better looking graphs -> AD_Project_PCA.ipynb

#--------------------------------------------------------------------------------------------------------
#Performs MDS analysis

#Transpose the PCA data matrix first
data_for_PCA2 <- as.data.frame(t(as.matrix(data_for_PCA_filtered)))

#Calcuate MDS
mds <- cmdscale(dist(data_for_PCA2))

#Samples representation of PCA
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.3) 

############################################################
#Creating Design Model Matrix
############################################################

#--------------------------------------------------------------------------------------------------------
#create the design model matrix. Here it's just by group, but if you have different batches, 
#you want to define batches and use ~0+group+batch
#design matrix is made to compare between different conditions
design <- model.matrix(~0+group, data = sampleTable)

#Rename the colnames
colnames(design) <- levels(group)

#Condition 1 is NDC, Condition 2 is V717I, Condition 3 is A79V, Condition 4 is N141I
#Here we compare the V717I mutation to the control, A79V mutation to the control
#N141I mutation to the control 
contr.matrix <- makeContrasts(V717IvsNDC = V717I - NDC,
                              A79VvsNDC = A79V - NDC,
                              N141IvsNDC = N141I - NDC,
                              levels = colnames(design))

levels= colnames(design)
print(levels)

#Load the contrast matrix
contr.matrix

#--------------------------------------------------------------------------------------------------------
#Lets look at an MDS plot of the data
par(mfrow=c(1,2))

#Color by group
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")

#View an MDS plot
plotMDS(filteredExpr2, labels = group, col = levels(col.group))

#Use Glimma to view an interactive MDS plot
glMDSPlot(filteredExpr2, labels=paste(group), groups=c(y$samples))

############################################################
#DEG Analysis by Limma Voom - eBayes Method
############################################################

#--------------------------------------------------------------------------------------------------------
#Run voom program to normalize the counts
#Two graphs show up: the left graph shows gene-wise square root residual
#Gene-wise square-root residual standard deviations are plotted against 
#average log-count. The right graph shows the different samples against the
#weights
#Run voomwithQualityWeights
v <- voomWithQualityWeights(y, design=design, plot=TRUE)

# fit the linear model
fit <- lmFit(v)
names(fit)

# apply your contrasts
cfit <- contrasts.fit(fit, contrasts=contr.matrix)

# eBayes method
efit <- eBayes(cfit)
plotSA(efit, main="Final model: Mean-variance trend")
dim(efit)

#See DEGs by eBayes method
summary(decideTests(efit))

# get topTable results for all 3 comparisons
res_limma_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)
res_limma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)
res_limma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

# save the limma-voom object for future analysis 
save(efit, file = "~/RNA-seq/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")
