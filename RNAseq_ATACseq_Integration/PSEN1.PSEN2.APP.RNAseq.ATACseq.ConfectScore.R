#RNA-seq/ATAC-seq confect ordering anylsis
#Valdes et. al 2023 Molecular Psychiatry Submission
#load packages
library("DiffBind")
library("ChIPseeker")
library("rtracklayer")
library("plyranges")
library("topconfects")
library("edgeR")
# load db
txdb.hg38ens <- loadDb("~/gmtfiles/txdb.hg38ens.sqlite")
# load PEREGRINE enhancers
load("~/gmtfiles/PEREGRINE.brain.enhancers.Rdata")
# load ATAC-seq data
load("~/RNA_ATAC/PSEN1.PSEN2.APP.ATACSeqDBAcount.Rdata")
load("~/ATAC-seq/ATACSeq.V717I.Rdata")
load("~/ATAC-seq/ATACSeq.A79V.Rdata")
load("~/ATAC-seq/ATACSeq.N141I.Rdata")
# load RNA-seq data
load("~/RNA-seq/PSEN12APP_GRCh38.99_txilsTPM.Rdata")
load("~/RNA-seq/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")
#
ATACSeq.V717I.full <- dba.report(ATACSeq.V717I, method = DBA_EDGER, th=1)
ATACSeq.A79V.full <- dba.report(ATACSeq.A79V, method = DBA_EDGER, th=1)
ATACSeq.N141I.full <- dba.report(ATACSeq.N141I, method = DBA_EDGER, th=1)
# Create a sample table
rna_dir <- "~/RNA-seq"
samples.RNA <- read.table(file.path(rna_dir, "PSEN1_PSEN2_APP_seqinfo.txt"), header = TRUE, stringsAsFactors=FALSE)
sampleTable <- data.frame(condition=factor(rep(c(samples.RNA$condition))))
# Create DGEList
y <- DGEList(txi_lsTPM$counts,
             lib.size = colSums(txi_lsTPM$counts),
             norm.factors = calcNormFactors(txi_lsTPM$counts),
             samples = samples.RNA$sample,
             group = samples.RNA$condition)
# Add the ENTREZID, SYMBOL for each ENSEMBL ID
library(org.Hs.eg.db)
#Create a Homo Sapiens annotation
Hs_ann <- AnnotationDbi::select(org.Hs.eg.db,
                                keys=rownames(y$counts),
                                columns=c("ENTREZID","SYMBOL"),
                                keytype="ENSEMBL",
                                multiVals="first")
# Remove duplicated terms
Hs_ann <- Hs_ann[!duplicated(Hs_ann[,1]),]
head(Hs_ann)
# Apply the annotation to your limma object "y"
y$genes <- Hs_ann
# View the library size for each sample
y$samples
# Filter lowly expressed genes via edgeR
keep = filterByExpr(y)
y <- y[keep,]
# Calculate normalization factors
y <- calcNormFactors(y, method = "TMM")
# Create the group and design - rename your conditions!
group <- factor(samples.RNA$condition,levels=c("NDC","V717I","A79V","N141I"))
# Define the model matrix 
design <- model.matrix(~0+group, data = sampleTable)
# Create the contrast matrix for DE
contr.matrix <- makeContrasts(V717IvNDC = groupV717I - groupNDC,
                              A79VvNDC = groupA79V - groupNDC,
                              N141IvNDC = groupN141I - groupNDC,
                              levels = design)
# See the contrast matrix
contr.matrix
# Can also run voomwithQualityWeights (see limma vignette) - unhash the next line if so
v <- voom(y, design=design, plot=TRUE)
# fit the linear model
fit <- lmFit(v, design)
# create contrast fit object for topconfects
cfit <- contrasts.fit(fit, contrasts = contr.matrix)

# Create DGEList for ATAC-seq using DiffBind and edgeR in order to derive confect scores from ATAC-seq data
samples <- ATACSeqDBAcount$samples
ATACcolData <- subset(samples, select = c(SampleID,Factor,Condition))
ATACcounts <- dba.peakset(ATACSeqDBAcount, bRetrieve=TRUE)
ATACcounts.df <- as.data.frame(ATACcounts)
ATACcountsonly <- subset(ATACcounts.df, select = -c(seqnames,start,end,width,strand))
ATAClimma <- SummarizedExperiment(assays = list(counts=ATACcountsonly),
                                  rowRanges = granges(ATACcounts),
                                  colData = ATACcolData)
ATAC.DGEList <- edgeR::SE2DGEList(ATAClimma)
# Create the group and design - rename your conditions!
condition.ATAC <- factor(ATACcolData$Condition)
# Create the design matrix
design.ATAC <- model.matrix(~0+condition.ATAC, data = colData(ATAClimma))
colnames(design.ATAC) <- levels(condition.ATAC)
# Normalize
ATAC.DGEList <- edgeR::calcNormFactors(ATAC.DGEList)

# V717I RNA-seq confects score
confects.RNA.V717I <- limma_confects(cfit,coef = "V717IvNDC",fdr = 0.05,step = 0.01,trend = FALSE,full = FALSE)
confects.RNA.V717I.filter <- confects.RNA.V717I$table
confects.RNA.V717I.filter <- confects.RNA.V717I.filter[,c(7,9,3)]
colnames(confects.RNA.V717I.filter)[3] <- "V717I.de_confect"
confects.RNA.V717I.filter$V717I.de_confect[is.na(confects.RNA.V717I.filter$V717I.de_confect)] <- 0
head(confects.RNA.V717I.filter)
# full RNA-seq statistics
RNASeq.V717I.full <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)
RNASeq.V717I.full <- RNASeq.V717I.full %>%
  dplyr::select(ENSEMBL, V717I.de_log2FC = de_log2FC, V717I.de_padj = de_padj)
dim(RNASeq.V717I.full)
# Remove duplicated ENSEMBL terms
RNASeq.V717I.full = subset(RNASeq.V717I.full, ENSEMBL != "" )
RNASeq.V717I.full = subset(RNASeq.V717I.full, ! duplicated(ENSEMBL))
dim(RNASeq.V717I.full)
dim(confects.RNA.V717I.filter)
head(RNASeq.V717I.full)
head(confects.RNA.V717I.filter)
rownames(confects.RNA.V717I.filter) <- confects.RNA.V717I.filter$ENSEMBL
RNASeq.V717I.full.confects <- merge(RNASeq.V717I.full,confects.RNA.V717I.filter)
head(RNASeq.V717I.full.confects)
rownames(RNASeq.V717I.full.confects) <- RNASeq.V717I.full.confects$ENSEMBL

# V717I ATAC-seq confects score
contr.matrix.ATAC.V717I <- makeContrasts(NDCvV717I = V717I - NDC,
                                         levels = design.ATAC)
# edgeR glm for ATAC-seq to run topconfects
glm.ATAC.V717I <- edgeR::glmQLFit(ATAC.DGEList, design.ATAC, dispersion = ATACSeq.V717I$contrasts[[1]]$edgeR$GLM$dispersion)
confects.ATAC.V717I <- edger_confects(glm.ATAC.V717I,contrast = contr.matrix.ATAC.V717I, fdr = 0.05,step = 0.01)
confects.ATAC.V717I.filter <- confects.ATAC.V717I$table[,c(7:11,6,3)]
colnames(confects.ATAC.V717I.filter)[c(6:7)] <- c("PeakID","V717I.da_confect")
head(confects.ATAC.V717I.filter)
confects.ATAC.V717I.filter[is.na(confects.ATAC.V717I.filter)] <- 0
confects.ATAC.V717I.ranges <- as_granges(confects.ATAC.V717I.filter)
names(confects.ATAC.V717I.ranges) <- confects.ATAC.V717I.ranges$PeakID
# Full ATAC-seq statistics 
ATACSeq.V717I.full <- dba.report(ATACSeq.V717I, method = DBA_EDGER, th=1)
ATACSeq.V717I.full$PeakID <- names(ATACSeq.V717I.full)
anno.V717I <- annotatePeak(ATACSeq.V717I.full, 
                           tssRegion=c(-1500,500), 
                           TxDb=txdb.hg38ens, 
                           level="gene", 
                           annoDb="org.Hs.eg.db",
                           overlap="TSS")
anno.V717I.ranges <- as_granges(anno.V717I)
# V717I
anno.V717I.ranges <- anno.V717I.ranges %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )
# Select key columns and convert back to granges
ATACpeaks.V717I <- as_granges(as.data.frame(anno.V717I.ranges)[,c(1:5,12,13,24:26,22)])
names(ATACpeaks.V717I) <- ATACpeaks.V717I$PeakID
ATACpeaks.V717I.confect <- ATACpeaks.V717I %>%
  plyranges::join_overlap_intersect(confects.ATAC.V717I.ranges)
ATACpeaks.V717I.confect <- ATACpeaks.V717I.confect %>%
  dplyr::select(PeakID = PeakID.x, annotation, V717I.da_log2FC = da_log2FC, V717I.da_padj = da_padj, V717I.da_confect, ENSEMBL, SYMBOL )
#
ATACpeaks.V717I.w <- ATACpeaks.V717I.confect[na.omit(match(rownames(RNASeq.V717I.full.confects), ATACpeaks.V717I.confect$ENSEMBL))]
# Add the differential expression statistics of genes ID'd by RNA-Seq and ATAC-Seq
mcols(ATACpeaks.V717I.w) <- RNASeq.V717I.full.confects[match(ATACpeaks.V717I.w$ENSEMBL, rownames(RNASeq.V717I.full.confects)),]
# View the intergrated results
ATACpeaks.V717I.w
# Overlap the dATAC with dE peak regions
ATACpeaks.V717I.final <- ATACpeaks.V717I.confect %>%
  join_overlap_left(ATACpeaks.V717I.w)

# Rename columns and filter
RNA.ATAC.V717I <- ATACpeaks.V717I.final %>%
  dplyr::mutate(ENSEMBL = ENSEMBL.x,
                SYMBOL = SYMBOL.x,
                V717I.da_log2FC,
                V717I.da_padj,
                V717I.da_confect,
                V717I.de_log2FC,
                V717I.de_padj,
                V717I.de_confect)
RNA.ATAC.V717I <- RNA.ATAC.V717I[,c(1,2,13,14,3:5,9:10,12)]
RNA.ATAC.V717I.sig <- RNA.ATAC.V717I %>%
  filter(V717I.de_padj < 0.05) %>%
  filter(abs(V717I.de_log2FC) > 0) %>%
  filter(V717I.da_padj < 0.05) %>%
  filter(abs(V717I.da_log2FC) > 0)
RNA.ATAC.V717I.sig.confects <- RNA.ATAC.V717I.sig %>%
  filter(abs(V717I.de_confect) > 0) %>%
  filter(abs(V717I.da_confect) > 0)
RNA.ATAC.V717I.sig.confects
#
names(RNA.ATAC.V717I) <- NULL
RNA.ATAC.V717I.df <- as.data.frame(unique(RNA.ATAC.V717I))
RNA.ATAC.V717I.df$V717I.de_confect[is.na(RNA.ATAC.V717I.df$V717I.de_confect)] <- 0
# calculate z-score
RNA.ATAC.V717I.df$V717I.zscore <-
  (RNA.ATAC.V717I.df$V717I.de_confect/stats::sd(
    RNA.ATAC.V717I.df$V717I.de_confect))*
  (RNA.ATAC.V717I.df$V717I.da_confect/
     stats::sd(RNA.ATAC.V717I.df$V717I.da_confect))
RNA.ATAC.V717I.zscore <- as_granges(RNA.ATAC.V717I.df)
names(RNA.ATAC.V717I.zscore) <- RNA.ATAC.V717I.zscore$PeakID
# final topconfects zscore
RNA.ATAC.V717I.zscore

# A79V RNA-seq confects score
confects.RNA.A79V <- limma_confects(cfit,coef = "A79VvNDC",fdr = 0.05,step = 0.01,trend = FALSE,full = FALSE)
confects.RNA.A79V.filter <- confects.RNA.A79V$table
confects.RNA.A79V.filter <- confects.RNA.A79V.filter[,c(7,9,3)]
colnames(confects.RNA.A79V.filter)[3] <- "A79V.de_confect"
confects.RNA.A79V.filter$A79V.de_confect[is.na(confects.RNA.A79V.filter$A79V.de_confect)] <- 0
head(confects.RNA.A79V.filter)
# full RNA-seq statistics
RNASeq.A79V.full <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)
RNASeq.A79V.full <- RNASeq.A79V.full %>%
  dplyr::select(ENSEMBL, A79V.de_log2FC = de_log2FC, A79V.de_padj = de_padj)
dim(RNASeq.A79V.full)
# Remove duplicated ENSEMBL terms
RNASeq.A79V.full = subset(RNASeq.A79V.full, ENSEMBL != "" )
RNASeq.A79V.full = subset(RNASeq.A79V.full, ! duplicated(ENSEMBL))
dim(RNASeq.A79V.full)
head(RNASeq.A79V.full)
RNASeq.A79V.full.confects <- merge(RNASeq.A79V.full,confects.RNA.A79V.filter)
head(RNASeq.A79V.full.confects)
rownames(RNASeq.A79V.full.confects) <- RNASeq.A79V.full.confects$ENSEMBL

# A79V ATAC-seq confects score
contr.matrix.ATAC.A79V <- makeContrasts(NDCvA79V = A79V - NDC,
                                        levels = design.ATAC)
glm.ATAC.A79V <- edgeR::glmQLFit(ATAC.DGEList, design.ATAC, dispersion = ATACSeq.A79V$contrasts[[1]]$edgeR$GLM$dispersion)
confects.ATAC.A79V <- edger_confects(glm.ATAC.A79V,contrast = contr.matrix.ATAC.A79V, fdr = 0.05,step = 0.01)
confects.ATAC.A79V.filter <- confects.ATAC.A79V$table[,c(7:11,6,3)]
colnames(confects.ATAC.A79V.filter)[c(6:7)] <- c("PeakID","A79V.da_confect")
head(confects.ATAC.A79V.filter)
confects.ATAC.A79V.filter[is.na(confects.ATAC.A79V.filter)] <- 0
confects.ATAC.A79V.ranges <- as_granges(confects.ATAC.A79V.filter)
names(confects.ATAC.A79V.ranges) <- confects.ATAC.A79V.ranges$PeakID
# Full ATAC-seq statistics
ATACSeq.A79V.full$PeakID <- names(ATACSeq.A79V.full)
anno.A79V <- annotatePeak(ATACSeq.A79V.full, 
                          tssRegion=c(-1500,500), 
                          TxDb=txdb.hg38ens, 
                          level="gene", 
                          annoDb="org.Hs.eg.db",
                          overlap="TSS")
anno.A79V.ranges <- as_granges(anno.A79V)
# A79V
anno.A79V.ranges <- anno.A79V.ranges %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )
# Select key columns and convert back to granges
ATACpeaks.A79V <- as_granges(as.data.frame(anno.A79V.ranges)[,c(1:5,12,13,24:26,22)])
names(ATACpeaks.A79V) <- ATACpeaks.A79V$PeakID
ATACpeaks.A79V.confect <- ATACpeaks.A79V %>%
  plyranges::join_overlap_intersect(confects.ATAC.A79V.ranges)
ATACpeaks.A79V.confect <- ATACpeaks.A79V.confect %>%
  dplyr::select(PeakID = PeakID.x, annotation, A79V.da_log2FC = da_log2FC, A79V.da_padj = da_padj, A79V.da_confect, ENSEMBL, SYMBOL )
#
ATACpeaks.A79V.w <- ATACpeaks.A79V.confect[na.omit(match(rownames(RNASeq.A79V.full.confects), ATACpeaks.A79V.confect$ENSEMBL))]
# Add the differential expression statistics of genes ID'd by RNA-Seq and ATAC-Seq
mcols(ATACpeaks.A79V.w) <- RNASeq.A79V.full.confects[match(ATACpeaks.A79V.w$ENSEMBL, rownames(RNASeq.A79V.full.confects)),]
# View the intergrated results
ATACpeaks.A79V.w
# Overlap the dATAC with dE peak regions
ATACpeaks.A79V.final <- ATACpeaks.A79V.confect %>%
  join_overlap_left(ATACpeaks.A79V.w)
RNA.ATAC.A79V <- ATACpeaks.A79V.final %>%
  dplyr::mutate(ENSEMBL = ENSEMBL.x,
                SYMBOL = SYMBOL.x,
                A79V.da_log2FC,
                A79V.da_padj,
                A79V.da_confect,
                A79V.de_log2FC,
                A79V.de_padj,
                A79V.de_confect)
RNA.ATAC.A79V <- RNA.ATAC.A79V[,c(1,2,13,14,3:5,9:10,12)]
RNA.ATAC.A79V.sig <- RNA.ATAC.A79V %>%
  filter(A79V.de_padj < 0.05) %>%
  filter(abs(A79V.de_log2FC) > 0) %>%
  filter(A79V.da_padj < 0.05) %>%
  filter(abs(A79V.da_log2FC) > 0)
RNA.ATAC.A79V.sig.confects <- RNA.ATAC.A79V.sig %>%
  filter(abs(A79V.de_confect) > 0) %>%
  filter(abs(A79V.da_confect) > 0)
RNA.ATAC.A79V.sig.confects
#
names(RNA.ATAC.A79V) <- NULL
RNA.ATAC.A79V.df <- as.data.frame(unique(RNA.ATAC.A79V))
RNA.ATAC.A79V.df$A79V.de_confect[is.na(RNA.ATAC.A79V.df$A79V.de_confect)] <- 0
RNA.ATAC.A79V.df$A79V.zscore <-
  (RNA.ATAC.A79V.df$A79V.de_confect/stats::sd(
    RNA.ATAC.A79V.df$A79V.de_confect))*
  (RNA.ATAC.A79V.df$A79V.da_confect/
     stats::sd(RNA.ATAC.A79V.df$A79V.da_confect))
RNA.ATAC.A79V.zscore <- as_granges(RNA.ATAC.A79V.df)
names(RNA.ATAC.A79V.zscore) <- RNA.ATAC.A79V.zscore$PeakID
# final topconfects zscore
RNA.ATAC.A79V.zscore

# N141I RNA-seq confect score
confects.RNA.N141I <- limma_confects(cfit,coef = "N141IvNDC",fdr = 0.05,step = 0.01,trend = FALSE,full = FALSE)
confects.RNA.N141I.filter <- confects.RNA.N141I$table
confects.RNA.N141I.filter <- confects.RNA.N141I.filter[,c(7,9,3)]
colnames(confects.RNA.N141I.filter)[3] <- "N141I.de_confect"
confects.RNA.N141I.filter$N141I.de_confect[is.na(confects.RNA.N141I.filter$N141I.de_confect)] <- 0
head(confects.RNA.N141I.filter)
# Full RNA-seq statistics
RNASeq.N141I.full <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)
RNASeq.N141I.full <- RNASeq.N141I.full %>%
  dplyr::select(ENSEMBL, N141I.de_log2FC = de_log2FC, N141I.de_padj = de_padj)
dim(RNASeq.N141I.full)
# Remove duplicated ENSEMBL terms
RNASeq.N141I.full = subset(RNASeq.N141I.full, ENSEMBL != "" )
RNASeq.N141I.full = subset(RNASeq.N141I.full, ! duplicated(ENSEMBL))
dim(RNASeq.N141I.full)
head(RNASeq.N141I.full)
RNASeq.N141I.full.confects <- merge(RNASeq.N141I.full,confects.RNA.N141I.filter)
head(RNASeq.N141I.full.confects)
rownames(RNASeq.N141I.full.confects) <- RNASeq.N141I.full.confects$ENSEMBL

# N141I ATAC-seq confect score
contr.matrix.ATAC.N141I <- makeContrasts(NDCvN141I = A141I - NDC,
                                         levels = design.ATAC)
glm.ATAC.N141I <- edgeR::glmQLFit(ATAC.DGEList, design.ATAC, dispersion = ATACSeq.N141I$contrasts[[1]]$edgeR$GLM$dispersion)
confects.ATAC.N141I <- edger_confects(glm.ATAC.N141I,contrast = contr.matrix.ATAC.N141I, fdr = 0.05,step = 0.01)
confects.ATAC.N141I.filter <- confects.ATAC.N141I$table[,c(7:11,6,3)]
colnames(confects.ATAC.N141I.filter)[c(6:7)] <- c("PeakID","N141I.da_confect")
head(confects.ATAC.N141I.filter)
confects.ATAC.N141I.filter[is.na(confects.ATAC.N141I.filter)] <- 0
confects.ATAC.N141I.ranges <- as_granges(confects.ATAC.N141I.filter)
names(confects.ATAC.N141I.ranges) <- confects.ATAC.N141I.ranges$PeakID

# ATAC-seq full statistics
ATACSeq.N141I.full$PeakID <- names(ATACSeq.N141I.full)
anno.N141I <- annotatePeak(ATACSeq.N141I.full, 
                           tssRegion=c(-1500,500), 
                           TxDb=txdb.hg38ens, 
                           level="gene", 
                           annoDb="org.Hs.eg.db",
                           overlap="TSS")
anno.N141I.ranges <- as_granges(anno.N141I)
# N141I
anno.N141I.ranges <- anno.N141I.ranges %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )
# Select key columns and convert back to granges
ATACpeaks.N141I <- as_granges(as.data.frame(anno.N141I.ranges)[,c(1:5,12,13,24:26,22)])
names(ATACpeaks.N141I) <- ATACpeaks.N141I$PeakID
ATACpeaks.N141I.confect <- ATACpeaks.N141I %>%
  plyranges::join_overlap_intersect(confects.ATAC.N141I.ranges)
ATACpeaks.N141I.confect <- ATACpeaks.N141I.confect %>%
  dplyr::select(PeakID = PeakID.x, annotation, N141I.da_log2FC = da_log2FC, N141I.da_padj = da_padj, N141I.da_confect, ENSEMBL, SYMBOL )
#
ATACpeaks.N141I.w <- ATACpeaks.N141I.confect[na.omit(match(rownames(RNASeq.N141I.full.confects), ATACpeaks.N141I.confect$ENSEMBL))]
# Add the differential expression statistics of genes ID'd by RNA-Seq and ATAC-Seq
mcols(ATACpeaks.N141I.w) <- RNASeq.N141I.full.confects[match(ATACpeaks.N141I.w$ENSEMBL, rownames(RNASeq.N141I.full.confects)),]
# View the intergrated results
ATACpeaks.N141I.w
# Overlap the dATAC with dE peak regions
ATACpeaks.N141I.final <- ATACpeaks.N141I.confect %>%
  join_overlap_left(ATACpeaks.N141I.w)
RNA.ATAC.N141I <- ATACpeaks.N141I.final %>%
  dplyr::mutate(ENSEMBL = ENSEMBL.x,
                SYMBOL = SYMBOL.x,
                N141I.da_log2FC,
                N141I.da_padj,
                N141I.da_confect,
                N141I.de_log2FC,
                N141I.de_padj,
                N141I.de_confect)
RNA.ATAC.N141I <- RNA.ATAC.N141I[,c(1,2,13,14,3:5,9:10,12)]
RNA.ATAC.N141I.sig <- RNA.ATAC.N141I %>%
  filter(N141I.de_padj < 0.05) %>%
  filter(abs(N141I.de_log2FC) > 0) %>%
  filter(N141I.da_padj < 0.05) %>%
  filter(abs(N141I.da_log2FC) > 0)
RNA.ATAC.N141I.sig.confects <- RNA.ATAC.N141I %>%
  filter(abs(N141I.de_confect) > 0) %>%
  filter(abs(N141I.da_confect) > 0)
RNA.ATAC.N141I.sig.confects
#
names(RNA.ATAC.N141I) <- NULL
RNA.ATAC.N141I.df <- as.data.frame(unique(RNA.ATAC.N141I))
RNA.ATAC.N141I.df$N141I.de_confect[is.na(RNA.ATAC.N141I.df$N141I.de_confect)] <- 0
RNA.ATAC.N141I.df$N141I.zscore <-
  (RNA.ATAC.N141I.df$N141I.de_confect/stats::sd(
    RNA.ATAC.N141I.df$N141I.de_confect))*
  (RNA.ATAC.N141I.df$N141I.da_confect/
     stats::sd(RNA.ATAC.N141I.df$N141I.da_confect))
RNA.ATAC.N141I.zscore <- as_granges(RNA.ATAC.N141I.df)
names(RNA.ATAC.N141I.zscore) <- RNA.ATAC.N141I.zscore$PeakID
# final topconfects zscore
RNA.ATAC.N141I.zscore

# Filter results
RNA.ATAC.V717I <- RNA.ATAC.V717I  %>%
  filter(abs(V717I.de_log2FC) > 0) %>%
  filter(abs(V717I.da_log2FC) > 0)
RNA.ATAC.A79V <- RNA.ATAC.A79V  %>%
  filter(abs(A79V.de_log2FC) > 0) %>%
  filter(abs(A79V.da_log2FC) > 0)
RNA.ATAC.N141I <- RNA.ATAC.N141I  %>%
  filter(abs(N141I.de_log2FC) > 0) %>%
  filter(abs(N141I.da_log2FC) > 0)
# Intersect all 3 mutations
RNA.ATAC.V717I.A79V <- RNA.ATAC.V717I.zscore %>%
  join_overlap_left(RNA.ATAC.A79V.zscore)
RNA.ATAC.all3 <- RNA.ATAC.V717I.A79V %>%
  join_overlap_left(RNA.ATAC.N141I.zscore)
# subset columns
RNA.ATAC.all3 <- RNA.ATAC.all3[,c(23:26,5:11,16:22,27:33)]
# Rename
RNA.ATAC.all3.zscore <- RNA.ATAC.all3
# Split into promoter, enhancer, and distal intergenic regions
RNA.ATAC.all3.zscore.promoter <- RNA.ATAC.all3.zscore[RNA.ATAC.all3.zscore$annotation == "Promoter"]
RNA.ATAC.all3.zscore.promoter
RNA.ATAC.all3.zscore.nonpromoter <- RNA.ATAC.all3.zscore[RNA.ATAC.all3.zscore$annotation != "Promoter"]
RNA.ATAC.all3.zscore.nonpromoter
# Define enhancer using PEREGRINE
RNA.ATAC.all3.zscore.enhancer <- RNA.ATAC.all3.zscore.nonpromoter %>%
  join_overlap_inner(PEREGRINE.brain.ranges)
RNA.ATAC.all3.zscore.enhancer <- unique(RNA.ATAC.all3.zscore.enhancer)
# Distal intergenic defined as non-promoter, non-(PEREGRINE) enhancer peaks
RNA.ATAC.all3.zscore.di.regions <- GenomicRanges::setdiff(RNA.ATAC.all3.zscore.nonpromoter,RNA.ATAC.all3.zscore.enhancer)
RNA.ATAC.all3.zscore.di <- subsetByOverlaps(RNA.ATAC.all3.zscore,RNA.ATAC.all3.zscore.di.regions)
RNA.ATAC.all3.zscore.di
# Find Significant peak/genes combos for all 3 mutations
# V717I
RNA.ATAC.all3.zscore.promoter.V717I.confect <- RNA.ATAC.all3.zscore.promoter %>%
  dplyr::filter(abs(V717I.zscore) > 0)
RNA.ATAC.all3.zscore.promoter.V717I.sig <- RNA.ATAC.all3.zscore.promoter %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05)
RNA.ATAC.all3.zscore.enhancer.V717I.confect <- RNA.ATAC.all3.zscore.enhancer %>%
  dplyr::filter(abs(V717I.zscore) > 0)
RNA.ATAC.all3.zscore.enhancer.V717I.sig <- RNA.ATAC.all3.zscore.enhancer %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05)
RNA.ATAC.all3.zscore.di.V717I.confect <- RNA.ATAC.all3.zscore.di %>%
  dplyr::filter(abs(V717I.zscore) > 0)
RNA.ATAC.all3.zscore.di.V717I.sig <- RNA.ATAC.all3.zscore.di %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05)
# A79V
RNA.ATAC.all3.zscore.promoter.A79V.confect <- RNA.ATAC.all3.zscore.promoter %>%
  dplyr::filter(abs(A79V.zscore) > 0)
RNA.ATAC.all3.zscore.promoter.A79V.sig <- RNA.ATAC.all3.zscore.promoter %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05)
RNA.ATAC.all3.zscore.enhancer.A79V.confect <- RNA.ATAC.all3.zscore.enhancer %>%
  dplyr::filter(abs(A79V.zscore) > 0)
RNA.ATAC.all3.zscore.enhancer.A79V.sig <- RNA.ATAC.all3.zscore.enhancer %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05)
RNA.ATAC.all3.zscore.di.A79V.confect <- RNA.ATAC.all3.zscore.di %>%
  dplyr::filter(abs(A79V.zscore) > 0)
RNA.ATAC.all3.zscore.di.A79V.sig <- RNA.ATAC.all3.zscore.di %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05)
# N141I
RNA.ATAC.all3.zscore.promoter.N141I.confect <- RNA.ATAC.all3.zscore.promoter %>%
  dplyr::filter(abs(N141I.zscore) > 0)
RNA.ATAC.all3.zscore.promoter.N141I.sig <- RNA.ATAC.all3.zscore.promoter %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
RNA.ATAC.all3.zscore.enhancer.N141I.confect <- RNA.ATAC.all3.zscore.enhancer %>%
  dplyr::filter(abs(N141I.zscore) > 0)
RNA.ATAC.all3.zscore.enhancer.N141I.sig <- RNA.ATAC.all3.zscore.enhancer %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
RNA.ATAC.all3.zscore.di.N141I.confect <- RNA.ATAC.all3.zscore.di %>%
  dplyr::filter(abs(N141I.zscore) > 0)
RNA.ATAC.all3.zscore.di.N141I.sig <- RNA.ATAC.all3.zscore.di %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
# re-merge back
promoter.V717I.A79V <- c(RNA.ATAC.all3.zscore.promoter.V717I.confect,RNA.ATAC.all3.zscore.promoter.A79V.confect)
promoter.all3 <- unique(c(promoter.V717I.A79V,RNA.ATAC.all3.zscore.promoter.N141I.confect))
enhancer.V717I.A79V <- c(RNA.ATAC.all3.zscore.enhancer.V717I.confect,RNA.ATAC.all3.zscore.enhancer.A79V.confect)
enhancer.all3 <- unique(c(enhancer.V717I.A79V,RNA.ATAC.all3.zscore.enhancer.N141I.confect))
di.V717I.A79V <- c(RNA.ATAC.all3.zscore.di.V717I.confect,RNA.ATAC.all3.zscore.di.A79V.confect)
di.all3 <- unique(c(di.V717I.A79V,RNA.ATAC.all3.zscore.di.N141I.confect))
# final
promoter.all3
enhancer.all3
di.all3
# export tables
write.table(promoter.all3, file = "PSEN1.PSEN2.APP.RNA.ATAC.confect.zscore.promoter.tsv",sep = "\t", row.names = FALSE)
write.table(enhancer.all3, file = "PSEN1.PSEN2.APP.RNA.ATAC.confect.zscore.enhancer.tsv",sep = "\t", row.names = FALSE)
write.table(di.all3, file = "PSEN1.PSEN2.APP.RNA.ATAC.confect.zscore.di.tsv",sep = "\t", row.names = FALSE)
#
# save results
RNA.ATAC.all3.zscore.V717I.confect <- RNA.ATAC.all3.zscore %>%
  dplyr::filter(abs(V717I.zscore) > 0)
RNA.ATAC.all3.zscore.A79V.confect <- RNA.ATAC.all3.zscore %>%
  dplyr::filter(abs(A79V.zscore) > 0)
RNA.ATAC.all3.zscore.N141I.confect <- RNA.ATAC.all3.zscore %>%
  dplyr::filter(abs(N141I.zscore) > 0)
write.table(RNA.ATAC.all3.zscore.V717I.confect, file = "RNA.ATAC.zscore.sigconfect.V717I.tsv",sep = "\t", row.names = FALSE)
write.table(RNA.ATAC.all3.zscore.A79V.confect, file = "RNA.ATAC.zscore.sigconfect.A79V.tsv",sep = "\t", row.names = FALSE)
write.table(RNA.ATAC.all3.zscore.N141I.confect, file = "RNA.ATAC.zscore.sigconfect.N141I.tsv",sep = "\t", row.names = FALSE)
save(RNA.ATAC.all3.zscore,file="PSEN1.PSEN2.APP.RNA.ATAC.zscore.Rdata")
write.table(RNA.ATAC.all3.zscore, file = "PSEN1.PSEN2.APP.RNA.ATAC.zscore.tsv",sep = "\t", row.names = FALSE)

# promoter
zscore.promoter.all3 <- zscore.promoter %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.promoter.all3
zscore.promoter.rest.regions <- setdiff(zscore.promoter,zscore.promoter.all3)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter,zscore.promoter.rest.regions)
#
zscore.promoter.APP.PSEN1 <- zscore.promoter.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05)
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.APP.PSEN1)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.APP.PSEN2 <- zscore.promoter.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.APP.PSEN2)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.PSEN1.PSEN2 <- zscore.promoter.rest %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.PSEN1.PSEN2)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.APP.only <- zscore.promoter.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) 
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.APP.only)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.PSEN1.only <- zscore.promoter.rest %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) 
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.PSEN1.only)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.PSEN2.only <- zscore.promoter.rest %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05) 
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.PSEN2.only)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
# by confects
zscore.promoter.all3 <- zscore.promoter %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.promoter.all3
zscore.promoter.rest.regions <- setdiff(zscore.promoter,zscore.promoter.all3)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter,zscore.promoter.rest.regions)
#
zscore.promoter.APP.PSEN1 <- zscore.promoter.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0)
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.APP.PSEN1)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.APP.PSEN2 <- zscore.promoter.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.APP.PSEN2)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.PSEN1.PSEN2 <- zscore.promoter.rest %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.PSEN1.PSEN2)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.APP.only <- zscore.promoter.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) 
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.APP.only)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.PSEN1.only <- zscore.promoter.rest %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) 
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.PSEN1.only)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.PSEN2.only <- zscore.promoter.rest %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0) 
zscore.promoter.rest.regions <- setdiff(zscore.promoter.rest,zscore.promoter.PSEN2.only)
zscore.promoter.rest <- subsetByOverlaps(zscore.promoter.rest,zscore.promoter.rest.regions)
#
zscore.promoter.all3 <- zscore.promoter.all3 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.promoter.APP.PSEN1 <- zscore.promoter.APP.PSEN1 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.promoter.APP.PSEN2 <- zscore.promoter.APP.PSEN2 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.promoter.PSEN1.PSEN2 <- zscore.promoter.PSEN1.PSEN2 %>%
  dplyr::arrange(desc(A79V.da_log2FC))
zscore.promoter.APP.only <- zscore.promoter.APP.only %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.promoter.PSEN1.only <- zscore.promoter.PSEN1.only %>%
  dplyr::arrange(desc(A79V.da_log2FC))
zscore.promoter.PSEN2.only <- zscore.promoter.PSEN2.only %>%
  dplyr::arrange(desc(N141I.da_log2FC))
# ordering
zscore.promoter.working <- as_granges((rbind(as.data.frame(zscore.promoter.all3),as.data.frame(zscore.promoter.APP.PSEN1))))
zscore.promoter.working <- as_granges((rbind(as.data.frame(zscore.promoter.working),as.data.frame(zscore.promoter.APP.PSEN2))))
zscore.promoter.working <- as_granges((rbind(as.data.frame(zscore.promoter.working),as.data.frame(zscore.promoter.PSEN1.PSEN2))))
zscore.promoter.working <- as_granges((rbind(as.data.frame(zscore.promoter.working),as.data.frame(zscore.promoter.APP.only))))
zscore.promoter.working <- as_granges((rbind(as.data.frame(zscore.promoter.working),as.data.frame(zscore.promoter.PSEN1.only))))
zscore.promoter.working <- as_granges((rbind(as.data.frame(zscore.promoter.working),as.data.frame(zscore.promoter.PSEN2.only))))
zscore.promoter.working <- zscore.promoter.working[!duplicated(zscore.promoter.working$PeakID),]
zscore.promoter.ordered <- zscore.promoter.working
# final ordered promoter RNA/ATAC confect z-score
zscore.promoter.ordered
#
# enhancer
zscore.enhancer.all3 <- zscore.enhancer %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.enhancer.all3
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer,zscore.enhancer.all3)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer,zscore.enhancer.rest.regions)
#
zscore.enhancer.APP.PSEN1 <- zscore.enhancer.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05)
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.APP.PSEN1)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.APP.PSEN2 <- zscore.enhancer.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.APP.PSEN2)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.PSEN1.PSEN2 <- zscore.enhancer.rest %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.PSEN1.PSEN2)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.APP.only <- zscore.enhancer.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) 
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.APP.only)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.PSEN1.only <- zscore.enhancer.rest %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) 
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.PSEN1.only)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.PSEN2.only <- zscore.enhancer.rest %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05) 
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.PSEN2.only)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
# by confects
zscore.enhancer.all3 <- zscore.enhancer %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.enhancer.all3
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer,zscore.enhancer.all3)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer,zscore.enhancer.rest.regions)
#
zscore.enhancer.APP.PSEN1 <- zscore.enhancer.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0)
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.APP.PSEN1)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.APP.PSEN2 <- zscore.enhancer.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.APP.PSEN2)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.PSEN1.PSEN2 <- zscore.enhancer.rest %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.PSEN1.PSEN2)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.APP.only <- zscore.enhancer.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) 
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.APP.only)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.PSEN1.only <- zscore.enhancer.rest %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) 
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.PSEN1.only)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.PSEN2.only <- zscore.enhancer.rest %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0) 
zscore.enhancer.rest.regions <- setdiff(zscore.enhancer.rest,zscore.enhancer.PSEN2.only)
zscore.enhancer.rest <- subsetByOverlaps(zscore.enhancer.rest,zscore.enhancer.rest.regions)
#
zscore.enhancer.all3 <- zscore.enhancer.all3 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.enhancer.APP.PSEN1 <- zscore.enhancer.APP.PSEN1 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.enhancer.APP.PSEN2 <- zscore.enhancer.APP.PSEN2 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.enhancer.PSEN1.PSEN2 <- zscore.enhancer.PSEN1.PSEN2 %>%
  dplyr::arrange(desc(A79V.da_log2FC))
zscore.enhancer.APP.only <- zscore.enhancer.APP.only %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.enhancer.PSEN1.only <- zscore.enhancer.PSEN1.only %>%
  dplyr::arrange(desc(A79V.da_log2FC))
zscore.enhancer.PSEN2.only <- zscore.enhancer.PSEN2.only %>%
  dplyr::arrange(desc(N141I.da_log2FC))
# ordering
zscore.enhancer.working <- as_granges((rbind(as.data.frame(zscore.enhancer.all3),as.data.frame(zscore.enhancer.APP.PSEN1))))
zscore.enhancer.working <- as_granges((rbind(as.data.frame(zscore.enhancer.working),as.data.frame(zscore.enhancer.APP.PSEN2))))
zscore.enhancer.working <- as_granges((rbind(as.data.frame(zscore.enhancer.working),as.data.frame(zscore.enhancer.PSEN1.PSEN2))))
zscore.enhancer.working <- as_granges((rbind(as.data.frame(zscore.enhancer.working),as.data.frame(zscore.enhancer.APP.only))))
zscore.enhancer.working <- as_granges((rbind(as.data.frame(zscore.enhancer.working),as.data.frame(zscore.enhancer.PSEN1.only))))
zscore.enhancer.working <- as_granges((rbind(as.data.frame(zscore.enhancer.working),as.data.frame(zscore.enhancer.PSEN2.only))))
zscore.enhancer.working <- zscore.enhancer.working[!duplicated(zscore.enhancer.working$PeakID),]
zscore.enhancer.ordered <- zscore.enhancer.working
# Final ordered enhancer RNA/ATAC confect zscore
zscore.enhancer.ordered
#
# di
zscore.di.all3 <- zscore.di %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.di.all3
zscore.di.rest.regions <- setdiff(zscore.di,zscore.di.all3)
zscore.di.rest <- subsetByOverlaps(zscore.di,zscore.di.rest.regions)
#
zscore.di.APP.PSEN1 <- zscore.di.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05)
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.APP.PSEN1)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.APP.PSEN2 <- zscore.di.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.APP.PSEN2)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.PSEN1.PSEN2 <- zscore.di.rest %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05)
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.PSEN1.PSEN2)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.APP.only <- zscore.di.rest %>%
  dplyr::filter(V717I.da_padj < 0.05) %>%
  dplyr::filter(V717I.de_padj < 0.05) 
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.APP.only)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.PSEN1.only <- zscore.di.rest %>%
  dplyr::filter(A79V.da_padj < 0.05) %>%
  dplyr::filter(A79V.de_padj < 0.05) 
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.PSEN1.only)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.PSEN2.only <- zscore.di.rest %>%
  dplyr::filter(N141I.da_padj < 0.05) %>%
  dplyr::filter(N141I.de_padj < 0.05) 
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.PSEN2.only)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
# by confects
zscore.di.all3 <- zscore.di %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.di.all3
zscore.di.rest.regions <- setdiff(zscore.di,zscore.di.all3)
zscore.di.rest <- subsetByOverlaps(zscore.di,zscore.di.rest.regions)
#
zscore.di.APP.PSEN1 <- zscore.di.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0)
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.APP.PSEN1)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.APP.PSEN2 <- zscore.di.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.APP.PSEN2)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.PSEN1.PSEN2 <- zscore.di.rest %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0)
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.PSEN1.PSEN2)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.APP.only <- zscore.di.rest %>%
  dplyr::filter(abs(V717I.da_confect) > 0) %>%
  dplyr::filter(abs(V717I.de_confect) > 0) 
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.APP.only)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.PSEN1.only <- zscore.di.rest %>%
  dplyr::filter(abs(A79V.da_confect) > 0) %>%
  dplyr::filter(abs(A79V.de_confect) > 0) 
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.PSEN1.only)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.PSEN2.only <- zscore.di.rest %>%
  dplyr::filter(abs(N141I.da_confect) > 0) %>%
  dplyr::filter(abs(N141I.de_confect) > 0) 
zscore.di.rest.regions <- setdiff(zscore.di.rest,zscore.di.PSEN2.only)
zscore.di.rest <- subsetByOverlaps(zscore.di.rest,zscore.di.rest.regions)
#
zscore.di.all3 <- zscore.di.all3 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.di.APP.PSEN1 <- zscore.di.APP.PSEN1 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.di.APP.PSEN2 <- zscore.di.APP.PSEN2 %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.di.PSEN1.PSEN2 <- zscore.di.PSEN1.PSEN2 %>%
  dplyr::arrange(desc(A79V.da_log2FC))
zscore.di.APP.only <- zscore.di.APP.only %>%
  dplyr::arrange(desc(V717I.da_log2FC))
zscore.di.PSEN1.only <- zscore.di.PSEN1.only %>%
  dplyr::arrange(desc(A79V.da_log2FC))
zscore.di.PSEN2.only <- zscore.di.PSEN2.only %>%
  dplyr::arrange(desc(N141I.da_log2FC))
# ordering
zscore.di.working <- as_granges((rbind(as.data.frame(zscore.di.all3),as.data.frame(zscore.di.APP.PSEN1))))
zscore.di.working <- as_granges((rbind(as.data.frame(zscore.di.working),as.data.frame(zscore.di.APP.PSEN2))))
zscore.di.working <- as_granges((rbind(as.data.frame(zscore.di.working),as.data.frame(zscore.di.PSEN1.PSEN2))))
zscore.di.working <- as_granges((rbind(as.data.frame(zscore.di.working),as.data.frame(zscore.di.APP.only))))
zscore.di.working <- as_granges((rbind(as.data.frame(zscore.di.working),as.data.frame(zscore.di.PSEN1.only))))
zscore.di.working <- as_granges((rbind(as.data.frame(zscore.di.working),as.data.frame(zscore.di.PSEN2.only))))
zscore.di.working <- zscore.di.working[!duplicated(zscore.di.working$PeakID),]
ç <- zscore.di.working
# Final ordered di RNA/ATAC confect zscore
# Final ordered enhancer RNA/ATAC confect zscore

#
# Final ordered RNA/ATAC zscored confect scores for heatmap 
zscore.promoter.ordered
zscore.enhancer.ordered
zscore.di.ordered
#save results as Rdata, tsv, and bed files
save(zscore.promoter.ordered,file="all3.zscore.promoter.ordered.Rdata")
save(zscore.enhancer.ordered,file="all3.zscore.enhancer.ordered.Rdata")
save(zscore.di.ordered,file="all3.zscore.di.ordered.Rdata")
write.table(zscore.promoter.ordered,file="all3.zscore.promoter.ordered.tsv",sep="\t",row.names = FALSE)
write.table(zscore.enhancer.ordered,file="all3.zscore.enhancer.ordered.tsv",sep="\t",row.names = FALSE)
write.table(zscore.di.ordered,file="all3.zscore.di.ordered.tsv",sep="\t",row.names = FALSE)
rtracklayer::export.bed(as.data.frame(zscore.promoter.ordered)[,c(1:3)],con = "all3zscore.promoter.ordered.bed",format = "bed")
rtracklayer::export.bed(as.data.frame(zscore.enhancer.ordered)[,c(1:3)],con = "all3zscore.enhancer.ordered.bed",format = "bed")
rtracklayer::export.bed(as.data.frame(zscore.di.ordered)[,c(1:3)],con = "all3zscore.di.ordered.bed",format = "bed")

