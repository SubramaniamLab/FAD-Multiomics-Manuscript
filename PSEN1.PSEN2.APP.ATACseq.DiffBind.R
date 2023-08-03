#ATAC-seq data processing pipeline
#Valdez et. al 2023 Molecular Psychiatry Submission

########################################################
#1. Differential binding for ATAC-seq data 
########################################################

#Install packages
BiocManager::install("DiffBind")
BiocManager::install("ChIPseeker")
BiocManager::install("rtracklayer")
BiocManager::install("plyranges")
BiocManager::install("GenomicFeatures")
BiocManager::install("devtools", force = TRUE) 
BiocManager::install("chipenrich")
#Load the installed packages 
library("DiffBind") 
library("ChIPseeker")
library("rtracklayer")
library("plyranges") 
library("dplyr")
library("ggplot2")
library("nVennR") 
library("GenomicFeatures")
library("data.table")  

#Read in the Samples 
#Set the base directory containing your files 
base_dir <- "~/ATAC-seq/"
base_dir <- getwd()
#Reading in the peak sets using the DiffBind format
samples = read.csv("~/ATAC-seq/ATAC_HMMRATAC.csv")

#Differential binding analysis (DBA) usually works with peak sets, which are sets of genomic intervals
#representing candidate protein binding sites
#Genomic interval consists of chromosome, start and end position, and score indicating confidence
#in, or strength of, the peak
#Create a differential binding analysis (DBA) object 
#Use DiffBind csv input format
ATACSeqDBA <- dba(sampleSheet="~/ATAC-seq/ATAC_HMMRATAC.csv")

#View the DBA object
ATACSeqDBA

#Plot a correlation heatmap for ATAC-Seq data 
#Using the data from the peak calls, a correlation heatmap can be generated which gives an
#initial clustering of the samples using the cross-correlations of each row of the binding matrix:
par(mar=c(1,1,1,1))   #get rid of 'Error in plot.new() : figure margins too large'
plot(ATACSeqDBA)
dev.off() #clear the plotting area 

plot(ATACSeqDBAconsensus)

#Calculate the overlap rate
ATACSeq.overlap.rate <- dba.overlap(ATACSeqDBA,mode=DBA_OLAP_RATE)
ATACSeq.overlap.rate

#View the masks (different conditions)
names(ATACSeqDBA$masks)

#Once the peaksets are read in, a merging function finds all overlapping peaks and derives a 
#single set of unique genomic intervals covering all the supplied peaks (a consensus peakset for the
#experiment).
#Find the peaks common (i.e. overlapping peaks) for each condition
#What are the peaks found between the consensus with minimum overlap of 3?
ATACSeqDBAcon <- dba.peakset(ATACSeqDBA,consensus = DBA_CONDITION, minOverlap = 3)
ATACSeqDBAcon

#Create a peaklist with common peaks within conditions (multiple replicates)
ATACConsensus <- ATACSeqDBAcon
ATACConsensus <-  dba(ATACConsensus, mask=ATACConsensus$masks$Consensus)

#Derivation of the consensus peakset -> represents the overall set of candidate binding sites to
#be used further in analysis
ConsensusPeaks <- dba.peakset(ATACConsensus, bRetrieve=TRUE)
ConsensusPeaks   #candidate binding sites

#Then Count the reads over common ATAC-Seq peak intervals
#DiffBind can use the supplied sequence read files to count how many reads overlap each interval 
#for each unique sample. By default, the peaks in the consensus peakset are re-centered and
#trimmed based on calculating their summits (point of greatest read overlap) in order
#to provide more standardized peak intervals.

#Perform differential expression on the peaks (i.e. obtain the reads)
#All the peaks possible -> consensus peaks
#score = TMM normalized (using edgeR) using ChIP read counts and effective library size
ATACSeqDBAcount <- dba.count(ATACSeqDBA,
                             peaks=ConsensusPeaks,
                             bParallel = FALSE,
                             score = DBA_SCORE_TMM_READS_EFFECTIVE,
                             summits = FALSE)

#View the chromatin peak intervals and frip (Fraction reads in peaks)
#Shows that all the samples are using the same, 53,598 length consensus peakset (i.e. chromatin
#peak interval).
#The second is labeled FRiP, which stands for Fraction of Reads in Peaks. 
#This is the proportion of reads for that sample that overlap a peak in the
#consensus peakset, and can be used to indicate which samples show more enrichment overall.
ATACSeqDBAcount

#Plot the correlation after factoring in reads within peaks
plot(ATACSeqDBAcount)
dev.off()

#Save as a Rdata file
save(ATACSeqDBAcount, file = "~/ATAC-seq/ATACSeqDBAcount.Rdata")

###################################
#2. Perform PCA for ATAC-seq data
###################################

#################################################################################
#2a. Get ATAC-seq PCA scores from Count Data
#################################################################################

#Perform differential expression on the peaks 

#Compute the PCA using dba.plotPCA function
pca <- dba.plotPCA(ATACSeqDBAcount,label=DBA_CONDITION, method = DBA_EDGER_BLOCK)

#Retreive the loading values from the PCA
pca_loadings <- pca$panel.args

#Loadings along the x-axis (PC1)
pca_loadings_x_PC1 <- data.table(pca_loadings[[1]]$x)

#Loadings along the x-axis (PC1)
pca_loadings_y_PC2 <- data.table(pca_loadings[[1]]$y)

#Bind the two data columns together
pca_loadings_combined <- data.table(cbind(pca_loadings_x_PC1, pca_loadings_y_PC2))

#Write the .csv file with all PCA scores 
write.csv(pca_loadings_combined, '~/ATAC-seq/ATAC-seq_PCA_scores_PSEN1_PSEN2_APP_NDC.csv', append =FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

#################################################################################
#2b. #Plot PCA of conditions of differential peaks
#################################################################################

pca <- dba.plotPCA(ATACSeqDBAcount,label=DBA_CONDITION, method = DBA_EDGER_BLOCK)
plot(pca)
dev.off()

###############################################################################
#3. Find DA regions using edgeR-DiffBind method and normalize for batch factor 
#APP-V717I
###############################################################################

###V717I Condition
ATACSeq.contrast.V717I <- dba.contrast(ATACSeqDBAcount,
                                       ATACSeqDBAcount$masks$V717I,
                                       ATACSeqDBAcount$masks$NDC,
                                       "V717I",
                                       "NDC")

#Performing the main differential analysis function for the V717 vs. NDC contrast
ATACSeq.V717I <- dba.analyze(ATACSeq.contrast.V717I,method = DBA_ALL_METHODS)

#Create report that is similar to topTable in limma-voom with differentially expressed peaks
#DBA_EDGER is the differential analysis method: edgeR 
ATACSeq.report.V717I <- dba.report(ATACSeq.V717I, method = DBA_EDGER)
ATACSeq.report.V717I

#Write report to file
ATACSeq.report.V717I.out <- as.data.frame(ATACSeq.report.V717I)
write.csv(ATACSeq.report.V717I.out, file="~/ATAC-seq/ATAC-seq_V717I_vs_NDC_DARs.csv", sep="\t", quote=F, row.names=T)

#Create report that is similar to topTable in limma-voom with all peaks
ATACSeq.report.V717I_full <- dba.report(ATACSeq.V717I, method = DBA_EDGER, th =1)
ATACSeq.report.V717I_full

#Write full report to file
ATACSeq.report.V717I.out_full <- as.data.frame(ATACSeq.report.V717I_full)
write.csv(ATACSeq.report.V717I.out_full, file="~/ATAC-seq/ATAC-seq_V717I_vs_NDC_DARs_full.csv", sep="\t", quote=F, row.names=T)
save(ATACSeq.report.V717I_full, file = "~/ATAC-seq/PSEN1.PSEN2.APP.ATACSeq.report.V717I.full.Rdata")

################################################################
### 3a. START HERE for HINT Formatting! (APP-V717I vs. NDC)
################################################################

#Make the APP-V717I report -> save as RData file -> for GimmeMotif Enrichment 
#16,652 differential peaks here
ATACSeq.report.V717I

#53,598 all possible total peaks here
ATACSeq.report.V717I_full

#save as Rdata 
save(ATACSeq.report.V717I, file = "~/ATAC-seq/ATACSeq.report.V717I.Rdata")


#Filter for adj. pval < 0.05
ATACSeq.report.V717I.filtered <- ATACSeq.report.V717I %>%
  dplyr::filter(FDR < 0.05)

#Convert to a data frame
ATACSeq.report.V717I.filtered_df <- as.data.frame(ATACSeq.report.V717I.filtered)
ATACSeq.report.V717I.out_full_df <- as.data.frame(ATACSeq.report.V717I.out_full)

#Select columns 1-3 (chr or seqnames, start, end)
ATACSeq.report.V717I.filtered_df <- ATACSeq.report.V717I.filtered_df[, c(1:3)]
ATACSeq.report.V717I.out_full_df <- ATACSeq.report.V717I.out_full[, c(1:3)]

#Export as a BED file, e.g. 210511.dATACpeaks_full.V717I.bed (to use for HINT!)
rtracklayer::export.bed(ATACSeq.report.V717I.out_full_df, con = "~/ATAC-seq/dATACpeaks_full.V717I.bed", format = "bed")


#Summary of results for each tool
contrasts.ATACSeq.V717I <- dba.show(ATACSeq.V717I, bContrasts=TRUE)
contrasts.ATACSeq.V717I

#Plot the V717I contrast
dev.off()
plot(ATACSeq.V717I, contrast= 1)

#Create a MA plot of V717I contrast with the EDGER method
#MA plots are a useful way to visualize the relationship between the overall binding level at
#each site and the magnitude of the change in binding enrichment between conditions, as
#well as the effect of normalization on data. Plot shows the V717I vs. NDC contrast
#16,814 differentially bound sites identified here
dba.plotMA(ATACSeq.V717I, method=DBA_EDGER)

#Create a Volcano plot of contrast 1 with the EDGER method
#Similar to MA plots, Volcano plots also highlight significantly differentially bound sites and
#show their fold changes. Here, however, the confidence statistic (FDR or p-value) is shown on
#a negative log scale, helping visualize the relationship between the magnitude of fold changes
#and the confidence that sites are differentially bound.
dba.plotVolcano(ATACSeq.V717I, method=DBA_EDGER)

#Create a pvalue plot of contrast 1 with the EDGER method (boxplots)
#Boxplots provide a way to view how read distributions differ between classes of binding sites
#where 16,814 differentially bound sites for the V717I vs. NDC contrast are identified.
pvals_V717I <- dba.plotBox(ATACSeq.V717I)

###############################################################################
#4. Find DA regions using edgeR-DiffBind method and normalize for batch factor 
#PSEN1-A79V
###############################################################################


###A79V Condition
ATACSeq.contrast.A79V <- dba.contrast(ATACSeqDBAcount,
                                       ATACSeqDBAcount$masks$A79V,
                                       ATACSeqDBAcount$masks$NDC,
                                       "A79V",
                                       "NDC")

#Performing the main differential analysis function for the A79V vs. NDC contrast
ATACSeq.A79V <- dba.analyze(ATACSeq.contrast.A79V,method = DBA_ALL_METHODS)

#Create report that is similar to topTable in limma-voom with differentially expressed peaks
#DBA_EDGER is the differential analysis method: edgeR 
#28772 differential peaks
ATACSeq.report.A79V <- dba.report(ATACSeq.A79V, method = DBA_EDGER)
ATACSeq.report.A79V

#Write report to file
ATACSeq.report.A79V.out <- as.data.frame(ATACSeq.report.A79V)
write.csv(ATACSeq.report.A79V.out, file="~/ATAC-seq/ATAC-seq_A79V_vs_NDC_DARs.csv", sep="\t", quote=F, row.names=T)

#Create report that is similar to topTable in limma-voom with all peaks
ATACSeq.report.A79V_full <- dba.report(ATACSeq.A79V, method = DBA_EDGER, th =1)
ATACSeq.report.A79V_full

#Write report to file
ATACSeq.report.A79V.out_full <- as.data.frame(ATACSeq.report.A79V_full)
write.csv(ATACSeq.report.A79V.out_full, file="~/ATAC-seq/ATAC-seq_A79V_vs_NDC_DARs_full.csv", sep="\t", quote=F, row.names=T)
save(ATACSeq.report.A79V_full, file = "~/ATAC-seq/PSEN1.PSEN2.APP.ATACSeq.report.A79V.full.Rdata")

##########################################################
#4a. START HERE for HINT Formatting! (PSEN1-A79V vs. NDC)
##########################################################

#Make the PSEN1-A79V report -> save as RData file -> for GimmeMotifs Enrichment
#28,772 differential peaks here
ATACSeq.report.A79V

#53,345 all possible peaks here
ATACSeq.report.A79V_full

#save as Rdata 
save(ATACSeq.report.A79V, file = "~/ATAC-seq/ATACSeq.report.A79V.Rdata")


#Filter for adj. pval < 0.05
ATACSeq.report.A79V.filtered <- ATACSeq.report.A79V %>%
  dplyr::filter(FDR < 0.05)

#Convert to a data frame
ATACSeq.report.A79V.filtered_df <- as.data.frame(ATACSeq.report.A79V.filtered)
ATACSeq.report.A79V.out_full_df <- as.data.frame(ATACSeq.report.A79V_full)

#Select columns 1-3 (chr or seqnames, start, end)
ATACSeq.report.A79V.filtered_df <- ATACSeq.report.A79V.filtered_df[, c(1:3)]
ATACSeq.report.A79V.out_full_df <- ATACSeq.report.A79V.out_full_df[, c(1:3)]


#Export as a BED file, e.g. dATACpeaks_full.A79V.bed (to use for HINT)
rtracklayer::export.bed(ATACSeq.report.A79V.out_full_df, con = "~/ATAC-seq/dATACpeaks_full.A79V.bed", format = "bed")


#Summary of results for each tool
contrasts.ATACSeq.A79V <- dba.show(ATACSeq.A79V, bContrasts=TRUE)
contrasts.ATACSeq.A79V

#Plot the A79V contrast
dev.off()
plot(ATACSeq.A79V, contrast= 1)

#Create a MA plot of A79V contrast with the EDGER method
#MA plots are a useful way to visualize the relationship between the overall binding level at
#each site and the magnitude of the change in binding enrichment between conditions, as
#well as the effect of normalization on data. Plot shows the A79V vs. NDC contrast
#28,772 differentially bound sites identified here
dba.plotMA(ATACSeq.A79V, method=DBA_EDGER)

#Create a Volcano plot of contrast 1 with the EDGER method
#Similar to MA plots, Volcano plots also highlight significantly differentially bound sites and
#show their fold changes. Here, however, the confidence statistic (FDR or p-value) is shown on
#a negative log scale, helping visualize the relationship between the magnitude of fold changes
#and the confidence that sites are differentially bound.
dba.plotVolcano(ATACSeq.A79V, method=DBA_EDGER)

#Create a pvalue plot of contrast 1 with the EDGER method (boxplots)
#Boxplots provide a way to view how read distributions differ between classes of binding sites
#where 16,814 differentially bound sites for the A79V vs. NDC contrast are identified.
pvals_A79V <- dba.plotBox(ATACSeq.A79V)

###############################################################################
#5. Find DA regions using edgeR-DiffBind method and normalize for batch factor 
#PSEN2-N141I
###############################################################################


###N141I Condition
ATACSeq.contrast.N141I <- dba.contrast(ATACSeqDBAcount,
                                      ATACSeqDBAcount$masks$N141I,
                                      ATACSeqDBAcount$masks$NDC,
                                      "N141I",
                                      "NDC")

#Performing the main differential analysis function for the N141I vs. NDC contrast 
ATACSeq.N141I <- dba.analyze(ATACSeq.contrast.N141I,method = DBA_ALL_METHODS)

#Create report that is similar to topTable in limma-voom with differentially expressed peaks
#DBA_EDGER is the differential analysis method: edgeR 
ATACSeq.report.N141I <- dba.report(ATACSeq.N141I, method = DBA_EDGER)
#19,721 peaks found
ATACSeq.report.N141I

#Write report to file
ATACSeq.report.N141I.out <- as.data.frame(ATACSeq.report.N141I)
write.csv(ATACSeq.report.N141I.out, file="~/ATAC-seq/ATAC-seq_N141I_vs_NDC_DARs.csv", sep="\t", quote=F, row.names=T)

#Create report that is similar to topTable in limma-voom with all peaks
ATACSeq.report.N141I_full <- dba.report(ATACSeq.N141I, method = DBA_EDGER, th =1)
ATACSeq.report.N141I_full

#Write report to file
ATACSeq.report.N141I.out_full <- as.data.frame(ATACSeq.report.N141I_full)
write.csv(ATACSeq.report.N141I.out_full, file="~/ATAC-seq_DARs/ATAC-seq_N141I_vs_NDC_DARs_full.csv", sep="\t", quote=F, row.names=T)
save(ATACSeq.report.N141I_full, file = "~/ATAC-seq/PSEN1.PSEN2.APP.ATACSeq.report.N141I.full.Rdata")

###########################################################
#5a. START HERE for HINT Formatting! (PSEN2-N141I vs. NDC)
###########################################################

#Make the PSEN2-N141I report -> save as Rdata file -> run GimmeMotifs motif enrichment
#19,721 peaks found
ATACSeq.report.N141I

#53,345 possible peaks found
ATACSeq.report.N141I_full

#save as Rdata
save(ATACSeq.report.N141I, file = "~/ATAC-seq/ATACSeq.report.N141I.Rdata")

#Filter for adj. pval < 0.05
ATACSeq.report.N141I.filtered <- ATACSeq.report.N141I %>%
  dplyr::filter(FDR < 0.05) 

#Convert to a data frame
ATACSeq.report.N141I.filtered_df <- as.data.frame(ATACSeq.report.N141I.filtered)
ATACSeq.report.N141I.out_full_df <- as.data.frame(ATACSeq.report.N141I_full)


#Select columns 1-3 (chr or seqnames, start, end)
ATACSeq.report.N141I.filtered_df <- ATACSeq.report.N141I.filtered_df[, c(1:3)]
ATACSeq.report.N141I.out_full_df <- ATACSeq.report.N141I.out_full_df[, c(1:3)]

#Export as a BED file, e.g. 210511.dATACpeaks_full.N141I.bed
rtracklayer::export.bed(ATACSeq.report.N141I.out_full_df, con = "~/ATAC-seq/dATACpeaks_full.N141I.bed", format = "bed")


#Summary of results for each tool
contrasts.ATACSeq.N141I <- dba.show(ATACSeq.N141I, bContrasts=TRUE)
contrasts.ATACSeq.N141I

#Plot the N141I contrast
dev.off()
plot(ATACSeq.N141I, contrast= 1)

#Create a MA plot of N141I contrast with the EDGER method
#MA plots are a useful way to visualize the relationship between the overall binding level at
#each site and the magnitude of the change in binding enrichment between conditions, as
#well as the effect of normalization on data. Plot shows the N141I vs. NDC contrast
#19,721 differentially bound sites identified here
dba.plotMA(ATACSeq.N141I, method=DBA_EDGER)

#Create a Volcano plot of contrast 1 with the EDGER method
#Similar to MA plots, Volcano plots also highlight significantly differentially bound sites and
#show their fold changes. Here, however, the confidence statistic (FDR or p-value) is shown on
#a negative log scale, helping visualize the relationship between the magnitude of fold changes
#and the confidence that sites are differentially bound.
dba.plotVolcano(ATACSeq.N141I, method=DBA_EDGER)

#Create a pvalue plot of contrast 1 with the EDGER method (boxplots)
#Boxplots provide a way to view how read distributions differ between classes of binding sites
#where 19,721 differentially bound sites for the A79V vs. NDC contrast are identified.
pvals_N141I <- dba.plotBox(ATACSeq.N141I)

#############################
#6. Obtain hg38 (human) genes
#############################

#========================================================================================================
#Perform wget from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz
#Then perform gunzip hg38.ensGene.gtf.gz
txdb.hg38ens <- makeTxDbFromGFF(file.path(base_dir, "hg38.ensGene.gtf"), format = "gtf")
saveDb(txdb.hg38ens,file = "txdb.hg38ens.sqlite")
txdb.hg38ens <- loadDb("txdb.hg38ens.sqlite")
txdb.hg38ens
keytypes(txdb.hg38ens)
seqlevels(txdb.hg38ens)


############################################################
##7. APP-V717I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION 
############################################################

#========================================================================================================
# Annotate the regions of differential ATAC accessibility - can vary tss region if desired

#da.padj = differential ATAC fdr adj pval
#de.padj = differential expression (RNA) fdr adj pval
anno.V717I <- annotatePeak(ATACSeq.report.V717I, 
                           tssRegion=c(-1000,500), 
                           TxDb=txdb.hg38ens, 
                           level="gene", 
                           annoDb="org.Hs.eg.db",
                           overlap="TSS")

#Have null values here
anno.V717I.ranges <- as_granges(anno.V717I)


#Filter based on adjusted p-value < 0.05
anno.V717I.filtered <- anno.V717I.ranges %>%
  dplyr::filter(FDR < 0.05) 

#Rename differential peaks as a different variable name
dATACpeaks.V717I <- anno.V717I.filtered

#save as Rdata
save(dATACpeaks.V717I, file = "~/ATAC-seq/dATACpeaks.V717I.annotated.Rdata")

# Convert to dataframe
anno.V717I <- as.data.frame(anno.V717I)


#Filter only for adjusted p-value
anno.V717I <- anno.V717I %>%
  filter(FDR < 0.05) 

library(dplyr)

# Rename key columns
anno.V717I <- anno.V717I %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )

# Rename key columns
anno.V717I.filtered <- anno.V717I.filtered %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )

############################################################
##7a. APP-V717I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION
#All Differential Peak Regions 
############################################################

# Select key columns and convert back to granges ##ERROR HERE 'subscript contains out-of-bounds indices'!!
dATACpeaks.V717I <- as_granges(anno.V717I.filtered[,c(7, 15:17, 18, 19, 20)])


# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.V717I) <- dATACpeaks.V717I$ENSEMBL

#All the peaks that are significant
#16,652 peaks are significant here
dATACpeaks.V717I

# Isolate the peaks with only ENSEMBL IDs as metadata for matching with gene expression data
#NEED TO CHANGE THIS
dATACpeaks.V717I.ranges <- as_granges(anno.V717I[,c(1:5,25)])

# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.V717I.ranges) <- dATACpeaks.V717I$ENSEMBL

# View the significant differential ATAC peaks
#16,652 differential peaks 
dATACpeaks.V717I.ranges

###########################################################################################
##7b. APP-V717I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION- All Differential Peak Regions 
#Promoter and Non-Promoter Regions
###########################################################################################

# Split into promoter regions and find only unique peak intervals
#5065 peaks in promoter
dATACpeaks.V717I.promoter <- dATACpeaks.V717I[dATACpeaks.V717I$annotation == "Promoter"]


# Split into non-promoter regions and find only unique peak intervals
#11,587 peaks in nonpromoter
dATACpeaks.V717I.nonpromoter <- dATACpeaks.V717I[dATACpeaks.V717I$annotation != "Promoter"]

write.csv(dATACpeaks.V717I.promoter, file="~/ATAC-seq/dATACpeaks.V717I.promoter.csv", sep="\t", quote=F, row.names=T)
write.csv(dATACpeaks.V717I.nonpromoter, file="~/ATAC-seq/dATACpeaks.V717I.nonpromoter.csv")

# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.V717I.promoter <- unique(dATACpeaks.V717I.promoter)

V717I.dATAC.promoter.up <- dATACpeaks.V717I.promoter %>%
  dplyr::filter(da_log2FC > 0)
V717I.dATAC.promoter.down <- dATACpeaks.V717I.promoter %>%
  dplyr::filter(da_log2FC < 0)

#Select the increased and decreased promoter accessibility for the APP-V717I mutation
#1226 peaks are differential in the enhancer (increased accessibility within the promoter region)
V717I.dATAC.promoter.up
#3839 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
V717I.dATAC.promoter.down

write.csv(dATACpeaks.V717I.promoter.up, file="~/ATAC-seq/dATACpeaks.V717I.promoter.up.csv", sep="\t", quote=F, row.names=T)
write.csv(dATACpeaks.V717I.promoter.down, file="~/ATAC-seq/dATACpeaks.V717I.promoter.down.csv", sep="\t", quote=F, row.names=T)


###########################################################################################
##7c. APP-V717I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION- All Differential Peak Regions 
#PEREGRINE Enhancer 
###########################################################################################


# Find PEREGINE enhancer associated peaks
enhancer.PEREGRINE <- as.data.frame(read.table("~/ATAC-seq/Enhancer_Data/PEREGRINE.brain.enhancer.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(enhancer.PEREGRINE)[1] <- "seqnames"
names(enhancer.PEREGRINE)[2] <- "start"
names(enhancer.PEREGRINE)[3] <- "end"

# Convert to 3 column bed
enhancer.PEREGRINE.ranges <- plyranges::as_granges(enhancer.PEREGRINE[,c(1:3)])

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
#17,848 peaks in the enhancer
dATACpeaks.V717I.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.V717I.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
#3174 peaks in the enhancer
dATACpeaks.V717I.enhancer.PEREGRINE <- unique(dATACpeaks.V717I.enhancer.PEREGRINE)
V717I.dATAC.enhancer.ranges.PEREGRINE <- plyranges::reduce_ranges(dATACpeaks.V717I.enhancer.PEREGRINE)
V717I.dATAC.enhancer.ranges.PEREGRINE.df <- as.data.frame(V717I.dATAC.enhancer.ranges.PEREGRINE)
V717I.dATAC.enhancer.ranges.PEREGRINE.df <- V717I.dATAC.enhancer.ranges.PEREGRINE.df[,c(1:3)]
rtracklayer::export.bed(V717I.dATAC.enhancer.ranges.PEREGRINE.df, con = "~/ATAC-seq/Enhancer_Data/V717I.dATAC.PEREGRINE.enhancer.inner.bed", format = "bed")

V717I.dATAC.enhancer.PEREGRINE.up <- dATACpeaks.V717I.enhancer.PEREGRINE %>%
  dplyr::filter(da_log2FC > 0)
V717I.dATAC.enhancer.PEREGRINE.down <- dATACpeaks.V717I.enhancer.PEREGRINE %>%
  dplyr::filter(da_log2FC < 0)

#Select the increased and decreased enhancer accessibility for the APP-V717I mutation
#1,927 peaks are differential in the enhancer (increased accesibility within the enhancer region)
V717I.dATAC.enhancer.PEREGRINE.up
#1,247 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
V717I.dATAC.enhancer.PEREGRINE.down

write.csv(V717I.dATAC.enhancer.PEREGRINE.up, file="~/ATAC-seq/PEREGRINE_Enhancer/dATAC.PEREGRINE.enhancer.up.V717I.csv")
write.csv(V717I.dATAC.enhancer.PEREGRINE.down, file="~/ATAC-seq/PEREGRINE_Enhancer/dATAC.PEREGRINE.enhancer.down.V717I.csv")


###########################################################################################
##8. APP-V717I Mutation - EBAYES METHOD WITH topTABLE (Integration with RNA-seq and ATAC-seq)
###########################################################################################

#========================================================================================================
# Annotate the regions of differential ATAC accessibility - can vary tss region if desired

#Before running the rest of the code, remove the dATACpeaks.V717I and dATACpeaks.V717I.ranges
#variables and rerun those two lines of code again

#Using the eBayes method with topTable for differential gene expression
sigGenes.V717I <- topTable(efit2, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

# Isolate the key attributes and rename the columns
sigGenes.V717I <- sigGenes.V717I %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.V717I)

# Remove duplicated ENSEMBL terms
sigGenes.V717I = subset(sigGenes.V717I, ENSEMBL != "" )
sigGenes.V717I = subset(sigGenes.V717I, ! duplicated(ENSEMBL))
dim(sigGenes.V717I)
head(sigGenes.V717I)

# Match the genes from gene expression data to the differentially accessible ATAC peaks
dATACpeaks.V717I.w <- dATACpeaks.V717I.ranges[na.omit(match(rownames(sigGenes.V717I), names(dATACpeaks.V717I.ranges)))]
# Add the differential expression statistics of genes ID's by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.V717I.w) <- sigGenes.V717I[match(names(dATACpeaks.V717I.w), rownames(sigGenes.V717I)),]

# Overlap the dATAC with dE peak regions
dATACpeaks.V717I <- dATACpeaks.V717I %>%
  join_overlap_left(dATACpeaks.V717I.w)

# Filter for both differential gene expression and differential accessibility
dATACpeaks.V717I.DEG  <- dATACpeaks.V717I  %>%
  dplyr::filter(de_padj < 0.05) %>%
  dplyr::filter(abs(de_log2FC) > 0) %>%
  dplyr::filter(da_padj < 0.05) %>%
  dplyr::filter(abs(da_log2FC) > 0)

# View the final granges object of dATAC peaks with differential gene expression
#761 differential peaks
dATACpeaks.V717I.DEG 

# Isolate the key statistics of each
dATACpeaks.V717I.DEG <- dATACpeaks.V717I.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.V717I.DEG
save(dATACpeaks.V717I.DEG, file ="~/ATAC-seq/dATACpeaks.V717I.DEG.Rdata")

###########################################################################################
##8a. APP-V717I Mutation - EBAYES METHOD WITH topTABLE (Integration with RNA-seq and ATAC-seq)
#Promoter and Non-promoter regions
###########################################################################################

# Split into promoter and non-promoter regions and find only unique peak intervals
dATACpeaks.V717I.DEG.promoter <- dATACpeaks.V717I.DEG[dATACpeaks.V717I.DEG$annotation == "Promoter"]
#287 differential peaks in the promoter region
dATACpeaks.V717I.DEG.promoter
dATACpeaks.V717I.DEG.nonpromoter <- dATACpeaks.V717I.DEG[dATACpeaks.V717I.DEG$annotation != "Promoter"]
#474 differential peaks in the nonpromoter region
dATACpeaks.V717I.DEG.nonpromoter

write.csv(dATACpeaks.V717I.DEG, file="~/ATAC-seq/dATACpeaks.V717I.DEG.csv")
write.csv(dATACpeaks.V717I.DEG.promoter, file="~/ATAC-seq/dATACpeaks.V717I.DEG.promoter.csv", row.names = TRUE)
write.csv(dATACpeaks.V717I.DEG.nonpromoter, file="~/ATAC-seq/dATACpeaks.V717I.DEG.nonpromoter.csv")

rtracklayer::export.bed(ATACSeq.report.V717I, con = "~/ATAC-seq/dATACpeaks.V717I.bed", format = "bed")

# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.V717I.DEG.promoter <- unique(dATACpeaks.V717I.DEG.promoter)

V717I.dATAC.promoter.DEG.up <- dATACpeaks.V717I.DEG.promoter %>%
  dplyr::filter(da_log2FC > 0) %>%
  dplyr::filter(de_log2FC > 0)
V717I.dATAC.promoter.DEG.down <- dATACpeaks.V717I.DEG.promoter %>%
  dplyr::filter(da_log2FC < 0) %>%
  dplyr::filter(de_log2FC < 0)

#Select the increased and decreased promoter accessibility for the APP-V717I mutation
#152 peaks are differential in the enhancer (increased accesibility within the promoter region)
V717I.dATAC.promoter.DEG.up
#100 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
V717I.dATAC.promoter.DEG.down

write.csv(V717I.dATAC.promoter.DEG.up, file="~/ATAC-seq/dATAC.promoter.DEG.up.V717I.csv")
write.csv(V717I.dATAC.promoter.DEG.down, file="~/ATAC-seq/dATAC.promoter.DEG.down.V717I.csv")

###########################################################################################
##8b. APP-V717I Mutation - EBAYES METHOD WITH topTABLE (Integration with RNA-seq and ATAC-seq)
#PEREGRINE enhancer regions 
###########################################################################################

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
dATACpeaks.V717I.DEG.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.V717I.DEG.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
dATACpeaks.V717I.DEG.enhancer.PEREGRINE <- unique(dATACpeaks.V717I.DEG.enhancer.PEREGRINE)
V717I.dATAC.enhancer.PEREGRINE.DEG.ranges <- plyranges::reduce_ranges(dATACpeaks.V717I.DEG.enhancer.PEREGRINE)
V717I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- as.data.frame(V717I.dATAC.TSS.enhancer.PEREGRINE.DEG.ranges)
V717I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- V717I.dATAC.TSS.enhancer.PEREGRINE.DEG.ranges.df[,c(1:3)]
rtracklayer::export.bed(V717I.dATAC.enhancer.PEREGRINE.DEG.ranges.df, con = "~/ATAC-seq/PEREGRINE_Enhancer/V717I.dATAC.enhancer.PEREGRINE.DEG.inner.bed", format = "bed")
V717I.dATAC.enhancer.PEREGRINE.DEG.up <- dATACpeaks.V717I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC > 0)
V717I.dATAC.enhancer.PEREGRINE.DEG.down <- dATACpeaks.V717I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC < 0)

write.csv(dATACpeaks.V717I.DEG.enhancer.PEREGRINE, file="~/ATAC-seq/PEREGRINE_Enhancer/dATAC.enhancer.PEREGRINE.DEG.V717I.csv")

#Select the increased and decreased enhancer accessibility for the APP-V717I mutation
#182 peaks are differential in the enhancer (increased accesibility within the enhancer region)
V717I.dATAC.enhancer.PEREGRINE.DEG.up
#92 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
V717I.dATAC.enhancer.PEREGRINE.DEG.down

write.csv(V717I.dATAC.enhancer.PEREGRINE.DEG.up, file="~/ATAC-seq/PEREGRINE_EnhancerdATAC.enhancer.PEREGRINE.DEG.up.V717I.csv")
write.csv(V717I.dATAC.enhancer.PEREGRINE.DEG.down, file="~/ATAC-seq/PEREGRINE_EnhancerdATAC.enhancer.PEREGRINE.DEG.down.V717I.csv")
write.csv(dATACpeaks.V717I.DEG.enhancer.PEREGRINE, file="~/ATAC-seq/PEREGRINE_EnhancerdATAC.enhancer.PEREGRINE.DEG.V717I.csv")


#############################################################
##9. PSEN1-A79V Mutation - WITHOUT RNA-SEQ DATA INTEGRATION
#############################################################

#========================================================================================================
# Annotate the regions of differential ATAC accessibility - can vary tss region if desired

#da.padj = differential ATAC fdr adj pval
#de.padj = differential expression (RNA) fdr adj pval
anno.A79V <- annotatePeak(ATACSeq.report.A79V, 
                           tssRegion=c(-1000,500), 
                           TxDb=txdb.hg38ens, 
                           level="gene", 
                           annoDb="org.Hs.eg.db",
                           overlap="TSS")

#Filter at FDR < 0.05
anno.A79V.ranges <- as_granges(anno.A79V)
anno.A79V.filtered <- anno.A79V.ranges %>%
  dplyr::filter(FDR < 0.05) 

#28,772 differential peaks
dATACpeaks.A79V <- anno.A79V.filtered

#save as Rdata
save(dATACpeaks.A79V, file = "~/ATAC-seq/dATACpeaks.A79V.Rdata")

# Convert to dataframe
anno.A79V <- as.data.frame(anno.A79V)


#Filtered those with FDR < 0.05
anno.A79V <- anno.A79V %>%
  filter(FDR < 0.05) 

# Rename key columns
anno.A79V <- anno.A79V %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )

# Rename key columns
anno.A79V.filtered <- anno.A79V.filtered %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )

####################################################################
##9a. PSEN1-A79V Mutation - WITHOUT RNA-SEQ DATA INTEGRATION 
#All Differential Peak Regions
####################################################################

# Select key columns and convert back to granges ##ERROR HERE 'subscript contains out-of-bounds indices'!!
#dATACpeaks.V717I <- as_granges(anno.V717I.filtered[,c(1:5,12,21,20,23,24,25)])
dATACpeaks.A79V <- as_granges(anno.A79V.filtered[,c(7, 15:17, 18, 19, 20)])

# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.A79V) <- dATACpeaks.A79V$ENSEMBL

#All the peaks that are significant
#28,772 differential peaks
dATACpeaks.A79V

# Isolate the peaks with only ENSEMBL IDs as metadata for matching with gene expression data
dATACpeaks.A79V.ranges <- as_granges(anno.A79V[,c(1:5,25)])

# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.A79V.ranges) <- dATACpeaks.A79V$ENSEMBL

# View the significant differential ATAC peaks
dATACpeaks.A79V.ranges

##########################################################################################
##9b. PSEN1-A79V Mutation - WITHOUT RNA-SEQ DATA INTEGRATION - All Differential Peak Regions
#Promoter and Non-promoter Regions
##########################################################################################

# Split into promoterregions and find only unique peak intervals
#7,595 differential peaks
dATACpeaks.A79V.promoter <- dATACpeaks.A79V[dATACpeaks.A79V$annotation == "Promoter"]

# Split into non-promoter regions and find only unique peak intervals
#21,177 differential peaks
dATACpeaks.A79V.nonpromoter <- dATACpeaks.A79V[dATACpeaks.A79V$annotation != "Promoter"]

write.csv(dATACpeaks.A79V, file="~/ATAC-seq/dATACpeaks.A79V.csv")
write.csv(dATACpeaks.A79V.promoter, file="~/ATAC-seq/dATACpeaks.A79V.promoter.csv", sep="\t", quote=F, row.names=T)

write.csv(dATACpeaks.A79V.nonpromoter, file="~/ATAC-seq/dATACpeaks.A79V.nonpromoter.csv")
save(ATACSeq.report.A79V, file = "~/ATAC-seq/ATACSeq.report.A79V.Rdata")

# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.A79V.promoter <- unique(dATACpeaks.A79V.promoter)

A79V.dATAC.promoter.up <- dATACpeaks.A79V.promoter %>%
  dplyr::filter(da_log2FC > 0)
A79V.dATAC.promoter.down <- dATACpeaks.A79V.promoter %>%
  dplyr::filter(da_log2FC < 0)

#Select the increased and decreased promoter accessibility for the PSEN1-A79V mutation
#2,435 peaks are differential in the enhancer (increased accesibility within the promoter region) (FDR < 0.05)
A79V.dATAC.promoter.up
#5,160 peaks are differential in the enhancer (decreased accessibility within the enchancer region) (FDR < 0.05)
A79V.dATAC.promoter.down

write.csv(A79V.dATAC.promoter.up, file="~/ATAC-seq/dATACpeaks.A79V.promoter.up.csv")
write.csv(A79V.dATAC.promoter.down, file="~/ATAC-seq/dATACpeaks.A79V.promoter.down.csv")

###########################################################################################
##9c. PSEN1-A79V Mutation - WITHOUT RNA-SEQ DATA INTEGRATION- All Differential Peak Regions 
#PEREGRINE Enhancer
###########################################################################################

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
#34.210 peaks in the enhancer
dATACpeaks.A79V.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.A79V.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
#3174 peaks in the enhancer
dATACpeaks.A79V.enhancer.PEREGRINE <- unique(dATACpeaks.A79V.enhancer.PEREGRINE)
A79V.dATAC.enhancer.ranges.PEREGRINE <- plyranges::reduce_ranges(dATACpeaks.A79V.enhancer.PEREGRINE)
A79V.dATAC.enhancer.ranges.PEREGRINE.df <- as.data.frame(A79V.dATAC.enhancer.ranges.PEREGRINE)
A79V.dATAC.enhancer.ranges.PEREGRINE.df <- A79V.dATAC.enhancer.ranges.PEREGRINE.df[,c(1:3)]
rtracklayer::export.bed(A79V.dATAC.enhancer.ranges.PEREGRINE.df, con = "~/ATAC-seq/Enhancer_Data/A79V.dATAC.PEREGRINE.enhancer.inner.bed", format = "bed")

A79V.dATAC.enhancer.PEREGRINE.up <- dATACpeaks.A79V.enhancer.PEREGRINE %>%
  dplyr::filter(da_log2FC > 0)
A79V.dATAC.enhancer.PEREGRINE.down <- dATACpeaks.A79V.enhancer.PEREGRINE %>%
  dplyr::filter(da_log2FC < 0)

#Select the increased and decreased enhancer accessibility for the APP-V717I mutation
#3,658 peaks are differential in the enhancer (increased accesibility within the enhancer region)
A79V.dATAC.enhancer.PEREGRINE.up
#2,236 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
A79V.dATAC.enhancer.PEREGRINE.down

write.csv(A79V.dATAC.enhancer.PEREGRINE.up, file="~/ATAC-seq/Enhancer_Data/dATAC.PEREGRINE.enhancer.up.A79V.csv")
write.csv(A79V.dATAC.enhancer.PEREGRINE.down, file="~/ATAC-seq/Enhancer_Data/dATAC.PEREGRINE.enhancer.down.A79V.csv")


##############################################################################################
##10. PSEN1-A79V Mutation - EBAYES METHOD WITH topTABLE (Integration of RNA-seq and ATAC-seq)
##############################################################################################

#========================================================================================================
# Annotate the regions of differential ATAC accessibility - can vary tss region if desired

#Before running the rest of the code, remove the dATACpeaks.V717I and dATACpeaks.V717I.ranges
#variables and rerun those two lines of code again

#Using the eBayes method with topTable for differential gene expression
sigGenes.A79V <- topTable(efit2, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

# Isolate the key attributes and rename the columns
sigGenes.A79V <- sigGenes.A79V %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.A79V)

# Remove duplicated ENSEMBL terms
sigGenes.A79V = subset(sigGenes.A79V, ENSEMBL != "" )
sigGenes.A79V = subset(sigGenes.A79V, ! duplicated(ENSEMBL))
dim(sigGenes.A79V)
head(sigGenes.A79V)

# Match the genes from gene expresssion data to the differentially accessible ATAC peaks
dATACpeaks.A79V.w <- dATACpeaks.A79V.ranges[na.omit(match(rownames(sigGenes.A79V), names(dATACpeaks.A79V.ranges)))]
# Add the differential expression statistics of genes ID's by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.A79V.w) <- sigGenes.A79V[match(names(dATACpeaks.A79V.w), rownames(sigGenes.A79V)),]

# Overlap the dATAC with dE peak regions
dATACpeaks.A79V <- dATACpeaks.A79V %>%
  join_overlap_left(dATACpeaks.A79V.w)


# Filter for both differential gene expression and differential accessibility
dATACpeaks.A79V.DEG  <- dATACpeaks.A79V  %>%
  dplyr::filter(de_padj < 0.05) %>%
  dplyr::filter(abs(de_log2FC) > 0) %>%
  dplyr::filter(da_padj < 0.05) %>%
  dplyr::filter(abs(da_log2FC) > 0)

# View the final granges object of dATAC peaks with differential gene expression
#1666 differential ATAC peaks
dATACpeaks.A79V.DEG

# Isolate the key statistics of each
dATACpeaks.A79V.DEG <- dATACpeaks.A79V.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.A79V.DEG
save(dATACpeaks.A79V.DEG, file ="210426.dATACpeaks.A79V.DEG.Rdata")

###############################################################################################
##10a. PSEN1-A79V Mutation - EBAYES METHOD WITH topTABLE (Integration of RNA-seq and ATAC-seq)
#Promoter and Non-promoter Regions
###############################################################################################

# Split into promoter regions and find only unique peak intervals
#627 differential peaks in promoter region
dATACpeaks.A79V.DEG.promoter <- dATACpeaks.A79V.DEG[dATACpeaks.A79V.DEG$annotation == "Promoter"]

# Split into non-promoter regions and find only unique peak intervals
#1039 differential peaks in nonpromoter region
dATACpeaks.A79V.DEG.nonpromoter <- dATACpeaks.A79V.DEG[dATACpeaks.A79V.DEG$annotation != "Promoter"]

write.csv(dATACpeaks.A79V.DEG, file="~/ATAC-seq/dATACpeaks.A79V.DEG.csv")
write.csv(dATACpeaks.A79V.DEG.promoter, file="~/ATAC-seq/dATACpeaks.A79V.DEG.promoter.csv", row.names = TRUE)
write.csv(dATACpeaks.A79V.DEG.nonpromoter, file="~/ATAC-seq/dATACpeaks.A79V.DEG.nonpromoter.csv")

rtracklayer::export.bed(ATACSeq.report.A79V, con = "~/ATAC-seq/dATACpeaks.A79V.bed", format = "bed")

# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.A79V.DEG.promoter <- unique(dATACpeaks.A79V.DEG.promoter)

A79V.dATAC.promoter.DEG.up <- dATACpeaks.A79V.DEG.promoter %>%
  dplyr::filter(da_log2FC > 0) %>%
  dplyr::filter(de_log2FC > 0)
A79V.dATAC.promoter.DEG.down <- dATACpeaks.A79V.DEG.promoter %>%
  dplyr::filter(da_log2FC < 0) %>%
  dplyr::filter(de_log2FC < 0)

#Select the increased and decreased promoter accessibility for the PSEN1-A79V mutation
#353 peaks are differential in the enhancer (increased accesibility within the promoter region)
A79V.dATAC.promoter.DEG.up
#173 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
A79V.dATAC.promoter.DEG.down

write.csv(A79V.dATAC.promoter.DEG.up, file="~/ATAC-seq/dATAC.promoter.DEG.up.A79V.csv")
write.csv(A79V.dATAC.promoter.DEG.down, file="~/ATAC-seq/dATAC.promoter.DEG.down.A79V.csv")


############################################################################################
##10b. PSEN1-A79V Mutation - EBAYES METHOD WITH topTABLE (Integration of RNA-seq and ATAC-seq)
#PEREGRINE enhancer regions
############################################################################################

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
dATACpeaks.A79V.DEG.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.A79V.DEG.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
dATACpeaks.A79V.DEG.enhancer.PEREGRINE <- unique(dATACpeaks.A79V.DEG.enhancer.PEREGRINE)
A79V.dATAC.enhancer.PEREGRINE.DEG.ranges <- plyranges::reduce_ranges(dATACpeaks.A79V.DEG.enhancer.PEREGRINE)
A79V.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- as.data.frame(A79V.dATAC.enhancer.PEREGRINE.DEG.ranges)
A79V.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- A79V.dATAC.enhancer.PEREGRINE.DEG.ranges.df[,c(1:3)]
rtracklayer::export.bed(A79V.dATAC.enhancer.PEREGRINE.DEG.ranges.df, con = "~/ATAC-seq/Enhancer_data/dATAC.enhancer.PEREGRINE.DEG.inner.bed", format = "bed")
A79V.dATAC.enhancer.PEREGRINE.DEG.up <- dATACpeaks.A79V.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC > 0)
A79V.dATAC.enhancer.PEREGRINE.DEG.down <- dATACpeaks.A79V.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC < 0)

write.csv(dATACpeaks.A79V.DEG.enhancer.PEREGRINE, file="~/ATAC-seq/Enhancer_data/dATAC.enhancer.PEREGRINE.DEG.A79V.csv")

#Select the increased and decreased enhancer accessibility for the PSEN1-A79V mutation
#422 peaks are differential in the enhancer (increased accesibility within the enhancer region)
A79V.dATAC.enhancer.PEREGRINE.DEG.up
#281 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
A79V.dATAC.enhancer.PEREGRINE.DEG.down

write.csv(A79V.dATAC.enhancer.PEREGRINE.DEG.up, file="~/ATAC-seq/Enhancer_data/dATAC.enhancer.PEREGRINE.DEG.up.A79V.csv")
write.csv(A79V.dATAC.enhancer.PEREGRINE.DEG.down, file="~/ATAC-seq/Enhancer_data/dATAC.enhancer.PEREGRINE.DEG.down.A79V.csv")


##############################################################
##11. PSEN2-N141I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION 
##############################################################

#========================================================================================================
# Annotate the regions of differential ATAC accessibility - can vary tss region if desired

#da.padj = differential ATAC fdr adj pval
#de.padj = differential expression (RNA) fdr adj pval
anno.N141I <- annotatePeak(ATACSeq.report.N141I, 
                          tssRegion=c(-1000,500), 
                          TxDb=txdb.hg38ens, 
                          level="gene", 
                          annoDb="org.Hs.eg.db",
                          overlap="TSS")


#Have null values here
anno.N141I.ranges <- as_granges(anno.N141I)


anno.N141I.filtered <- anno.N141I.ranges %>%
  dplyr::filter(FDR < 0.05) 

#19,721 peaks found here
dATACpeaks.N141I <- anno.N141I.filtered

#save as Rdata
save(dATACpeaks.N141I, file = "~/ATAC-seq/dATACpeaks.N141I.Rdata")

# Convert to dataframe
anno.N141I <- as.data.frame(anno.N141I)

#Filter based on FDR < 0.05
anno.N141I <- anno.N141I %>%
  filter(FDR < 0.05) 

# Rename key columns
anno.N141I <- anno.N141I %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )

# Rename key columns
anno.N141I.filtered <- anno.N141I.filtered %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )

##########################################################################################
##11a. PSEN2-N141I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION - All Differential Peak Regions
##########################################################################################


# Select key columns and convert back to granges ##ERROR HERE 'subscript contains out-of-bounds indices'!!
dATACpeaks.N141I <- as_granges(anno.N141I.filtered[,c(7, 15:17, 18, 19, 20)])

# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.N141I) <- dATACpeaks.N141I$ENSEMBL

#All the peaks that are significant
#19,721 differential peaks found
dATACpeaks.N141I

# Isolate the peaks with only ENSEMBL IDs as metadata for matching with gene expression data
dATACpeaks.N141I.ranges <- as_granges(anno.N141I[,c(1:5,25)])

# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.N141I.ranges) <- dATACpeaks.N141I$ENSEMBL

# View the significant differential ATAC peaks
dATACpeaks.N141I.ranges

##########################################################################################
##11b. PSEN2-N141I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION - All Differential Peak Regions
#Promoter and Non-promoter Regions
##########################################################################################

# Split into promoter and non-promoter regions and find only unique peak intervals
#5,543 differential peaks found
dATACpeaks.N141I.promoter <- dATACpeaks.N141I[dATACpeaks.N141I$annotation == "Promoter"]

# Split into promoter and non-promoter regions and find only unique peak intervals
#14,178 differential peaks found in nonpromoter
dATACpeaks.N141I.nonpromoter <- dATACpeaks.N141I[dATACpeaks.N141I$annotation != "Promoter"]

#Save individual files
write.csv(dATACpeaks.N141I, file="~/ATAC-seq/dATACpeaks.N141I.csv")
write.csv(dATACpeaks.N141I.promoter, file="~/ATAC-seq/dATACpeaks.N141I.promoter.csv", sep="\t", quote=F, row.names=T)
write.csv(dATACpeaks.N141I.nonpromoter, file="~/ATAC-seq/dATACpeaks.N141I.nonpromoter.csv")
save(ATACSeq.report.N141I, file = "~/ATAC-seq/ATACSeq.report.N141I.Rdata")

# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.N141I.promoter <- unique(dATACpeaks.A141I.promoter)

N141I.dATAC.promoter.up <- dATACpeaks.N141I.promoter %>%
  dplyr::filter(da_log2FC > 0)
N141I.dATAC.promoter.down <- dATACpeaks.N141I.promoter %>%
  dplyr::filter(da_log2FC < 0)

#Select the increased and decreased promoter accessibility for the PSEN1-A79V mutation
#1547 peaks are differential in the enhancer (increased accessibility within the promoter region) (FDR < 0.05)
N141I.dATAC.promoter.up
#3996 peaks are differential in the enhancer (decreased accessibility within the enhancer reigon) (FDR < 0.05)
N141I.dATAC.promoter.down

write.csv(N141I.dATAC.promoter.up, file="~/ATAC-seq/dATACpeaks.N141I.promoter.up.csv")
write.csv(N141I.dATAC.promoter.down, file="~/ATAC-seq/dATACpeaks.N141I.promoter.down.csv")


###########################################################################################
##11c. PSEN2-N141I Mutation - WITHOUT RNA-SEQ DATA INTEGRATION- All Differential Peak Regions 
#PEREGRINE Enhancer
###########################################################################################

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
#21,966 peaks in the enhancer
dATACpeaks.N141I.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.N141I.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
#3174 peaks in the enhancer
dATACpeaks.N141I.enhancer.PEREGRINE <- unique(dATACpeaks.N141I.enhancer.PEREGRINE)
N141I.dATAC.enhancer.ranges.PEREGRINE <- plyranges::reduce_ranges(dATACpeaks.N141I.enhancer.PEREGRINE)
N141I.dATAC.enhancer.ranges.PEREGRINE.df <- as.data.frame(N141I.dATAC.enhancer.ranges.PEREGRINE)
N141I.dATAC.enhancer.ranges.PEREGRINE.df <- N141I.dATAC.enhancer.ranges.PEREGRINE.df[,c(1:3)]
rtracklayer::export.bed(N141I.dATAC.enhancer.ranges.PEREGRINE.df, con = "~/ATAC-seq/Enhancer_Data/N141I.dATAC.PEREGRINE.enhancer.inner.bed", format = "bed")

N141I.dATAC.enhancer.PEREGRINE.up <- dATACpeaks.A141I.enhancer.PEREGRINE %>%
  dplyr::filter(da_log2FC > 0)
N141I.dATAC.enhancer.PEREGRINE.down <- dATACpeaks.A141I.enhancer.PEREGRINE %>%
  dplyr::filter(da_log2FC < 0)

#Select the increased and decreased enhancer accessibility for the APP-V717I mutation
#2,385 peaks are differential in the enhancer (increased accesibility within the enhancer region)
N141I.dATAC.enhancer.PEREGRINE.up
#1,506 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
N141I.dATAC.enhancer.PEREGRINE.down

write.csv(N141I.dATAC.enhancer.PEREGRINE.up, file="~/ATAC-seq/Enhancer_Data/dATAC.PEREGRINE.enhancer.up.N141I.csv")
write.csv(N141I.dATAC.enhancer.PEREGRINE.down, file="~/ATAC-seq/Enhancer_Data/dATAC.PEREGRINE.enhancer.down.N141I.csv")


#############################################################################################
##12. PSEN2-N141I Mutation - EBAYES METHOD WITH topTABLE (Integration of RNA-seq and ATAC-seq)
#############################################################################################

#========================================================================================================
# Annotate the regions of differential ATAC accessibility - can vary tss region if desired


#Before running the rest of the code, remove the dATACpeaks.V717I and dATACpeaks.V717I.ranges
#variables and rerun those two lines of code again

#Using the eBayes method with topTable for differential gene expression
sigGenes.N141I <- topTable(efit2, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

# Isolate the key attributes and rename the columns
sigGenes.N141I <- sigGenes.N141I %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.N141I)

# Remove duplicated ENSEMBL terms
sigGenes.N141I = subset(sigGenes.N141I, ENSEMBL != "" )
sigGenes.N141I = subset(sigGenes.N141I, ! duplicated(ENSEMBL))
dim(sigGenes.N141I)
head(sigGenes.N141I)

#############################################################################################
##12a. PSEN2-N141I Mutation - EBAYES METHOD WITH topTABLE (Integration of RNA-seq and ATAC-seq)
#All Differential peak regions
#############################################################################################

# Match the genes from gene expresssion data to the differentially accessible ATAC peaks
dATACpeaks.N141I.w <- dATACpeaks.N141I.ranges[na.omit(match(rownames(sigGenes.N141I), names(dATACpeaks.N141I.ranges)))]
# Add the differential expression statistics of genes ID's by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.N141I.w) <- sigGenes.N141I[match(names(dATACpeaks.N141I.w), rownames(sigGenes.N141I)),]

# Overlap the dATAC with dE peak regions
dATACpeaks.N141I <- dATACpeaks.N141I %>%
  join_overlap_left(dATACpeaks.N141I.w)


# Filter for both differential gene expression and differential accessibility
dATACpeaks.N141I.DEG  <- dATACpeaks.N141I  %>%
  dplyr::filter(de_padj < 0.05) %>%
  dplyr::filter(abs(de_log2FC) > 0) %>%
  dplyr::filter(da_padj < 0.05) %>%
  dplyr::filter(abs(da_log2FC) > 0)

# View the final granges object of dATAC peaks with differential gene expression
#898 differential peaks with differential gene expression
dATACpeaks.N141I.DEG

# Isolate the key statistics of each
dATACpeaks.N141I.DEG <- dATACpeaks.N141I.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.N141I.DEG
save(dATACpeaks.N141I.DEG, file ="210426.dATACpeaks.N141I.DEG.Rdata")

#############################################################################################
##12b. PSEN2-N141I Mutation - EBAYES METHOD WITH topTABLE (Integration of RNA-seq and ATAC-seq)
#Promoter and Non-promoter regions
#############################################################################################


# Split into promoter and non-promoter regions and find only unique peak intervals
dATACpeaks.N141I.DEG.promoter <- dATACpeaks.N141I.DEG[dATACpeaks.N141I.DEG$annotation == "Promoter"]
#259 differential peaks in promoter region
dATACpeaks.N141I.DEG.promoter
dATACpeaks.N141I.DEG.nonpromoter <- dATACpeaks.N141I.DEG[dATACpeaks.N141I.DEG$annotation != "Promoter"]
#639 differential peaks in nonpromoter region
dATACpeaks.N141I.DEG.nonpromoter

write.csv(dATACpeaks.N141I.DEG, file="~/ATAC-seq/dATACpeaks.N141I.DEG.csv")
write.csv(dATACpeaks.N141I.DEG.promoter, file="~/ATAC-seq/dATACpeaks.N141I.DEG.promoter.csv", row.names = TRUE)
write.csv(dATACpeaks.N141I.DEG.nonpromoter, file="~/ATAC-seq/dATACpeaks.N141I.DEG.nonpromoter.csv")

rtracklayer::export.bed(ATACSeq.report.N141I, con = "~/ATAC-seq/dATACpeaks.N141I.bed", format = "bed")

# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.N141I.DEG.promoter <- unique(dATACpeaks.N141I.DEG.promoter)


N141I.dATAC.promoter.DEG.up <- dATACpeaks.N141I.DEG.promoter %>%
  dplyr::filter(da_log2FC > 0) %>%
  dplyr::filter(de_log2FC > 0)
N141I.dATAC.promoter.DEG.down <- dATACpeaks.N141I.DEG.promoter %>%
  dplyr::filter(da_log2FC < 0) %>%
  dplyr::filter(de_log2FC < 0)

#Select the increased and decreased promoter accessibility for the PSEN1-A79V mutation
#133 peaks are differential in the enhancer (increased accesibility within the promoter region)
N141I.dATAC.promoter.DEG.up
#86 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
N141I.dATAC.promoter.DEG.down

write.csv(N141I.dATAC.promoter.DEG.up, file="~/ATAC-seq/dATAC.promoter.DEG.up.N141I.csv")
write.csv(N141I.dATAC.promoter.DEG.down, file="~/ATAC-seq/dATAC.promoter.DEG.down.N141I.csv")


#############################################################################################
##12c. PSEN2-N141I Mutation - EBAYES METHOD WITH topTABLE (Integration of RNA-seq and ATAC-seq)
#PEREGRINE Enhancer Regions
#############################################################################################

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
dATACpeaks.N141I.DEG.enhancer <- plyranges::join_overlap_inner(dATACpeaks.N141I.DEG.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
dATACpeaks.N141I.DEG.enhancer.PEREGRINE <- unique(dATACpeaks.N141I.DEG.enhancer.PEREGRINE)
N141I.dATAC.enhancer.PEREGRINE.DEG.ranges <- plyranges::reduce_ranges(dATACpeaks.N141I.DEG.enhancer.PEREGRINE)
N141I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- as.data.frame(N141I.dATAC.enhancer.PEREGRINE.DEG.ranges)
N141I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- N141I.dATAC.enhancer.PEREGRINE.DEG.ranges.df[,c(1:3)]
rtracklayer::export.bed(N141I.dATAC.enhancer.PEREGRINE.DEG.ranges.df, con = "~/ATAC-seq/Enhancer_data/N141I.dATAC.enhancer.PEREGRINE.DEG.inner.bed", format = "bed")
N141I.dATAC.enhancer.PEREGRINE.DEG.up <- dATACpeaks.N141I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC > 0)
N141I.dATAC.enhancer.PEREGRINE.DEG.down <- dATACpeaks.N141I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC < 0)

write.csv(dATACpeaks.N141I.DEG.enhancer, file="~/ATAC-seq/Enhancer_data/dATAC.enhancer.PEREGRINE.DEG.N141I.csv")

#Select the increased and decreased enhancer accessibility for the PSEN1-A79V mutation
#249 peaks are differential in the enhancer (increased accesibility within the enhancer region)
N141I.dATAC.enhancer.PEREGRINE.DEG.up
#142 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
N141I.dATAC.enhancer.PEREGRINE.DEG.down

write.csv(N141I.dATAC.enhancer.PEREGRINE.DEG.up, file="~/ATAC-seq/Enhancer_data/dATAC.enhancer.PEREGRINE.DEG.up.N141I.csv")
write.csv(N141I.dATAC.enhancer.PEREGRINE.DEG.down, file="~/ATAC-seq/Enhancer_data/dATAC.enhancer.PEREGRINE.DEG.down.N141I.csv")

