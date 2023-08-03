#ATAC-seq GimmeMotifs Maelstrom preparation pipeline
#Valdes et. al 2023 Molecular Psychiatry Submission
#load packages
library("DiffBind")
library("ChIPseeker")
library("rtracklayer")
library("plyranges")
#load hg38
txdb.hg38ens <- loadDb("~/gmtfiles/txdb.hg38ens.sqlite")
#load PEREGRINE enhancers
load("~/gmtfiles/PEREGRINE.brain.enhancers.Rdata")
# ATACSeqDBA simple samplesheet to create output for GimmeMotifs
ATACSeqDBA <- dba(sampleSheet="~/PSEN1_PSEN2_APP/ATAC-seq/ATACSeqDBA.summit.samples.3col.csv")
ATACSeqDBA$config$cores
ATACSeqDBA$config$cores <- 8
ATACSeqDBAcount500 <- dba.count(ATACSeqDBA, summits=250,bParallel=FALSE)
save(ATACSeqDBAcount500, file = "PSEN1.PSEN2.APP.ATAC.summits.dbacount500.Rdata")
# View the chromatin peak intervals and frip (Fraction reads in peaks)
ATACSeqDBAcount <- ATACSeqDBAcount500
# Plot the correlation after factoring in reads within peaks
plot(ATACSeqDBAcount)
# Retrieve the TMM normalized count matrix
ATACcounts <- dba.peakset(ATACSeqDBAcount, bRetrieve=TRUE)
ATACcounts.df <- as.data.frame(ATACcounts)
ATACcountsonly <- subset(ATACcounts.df, select = -c(seqnames,start,end,width,strand))
# Create Summarized Experiment of ATAC-seq data
ATAClimma <- SummarizedExperiment(assays = list(ATACcountsonly),
                                  rowRanges = granges(ATACcounts),
                                  colData = ATACSeqDBAcount$samples)
# Create the group and design - rename your conditions!
condition <- factor(ATACSeqDBAcount$Condition)
# Create the design matrix
design <- model.matrix(~0+Condition, data = colData(ATAClimma))
# Create the contrast matrix
contr.matrix <- makeContrasts(NDCvV717I = ConditionV717I - ConditionNDC,
                              NDCvA79V = ConditionA79V - ConditionNDC,
                              NDCvN141I = ConditionA141I - ConditionNDC,
                              levels = design)
# run voom to log2 normalize 
v <- voom(assay(ATAClimma), design)
# Find consensus peakset
ATACSeqDBAcon <- dba.peakset(ATACSeqDBA,consensus = DBA_CONDITION, minOverlap = 3)
ATACSeqDBAcon
# Create a peaklist with common peaks within conditions (multiple replicates)
ATACConsensus <- ATACSeqDBAcon
ATACConsensus <-  dba(ATACConsensus, mask=ATACConsensus$masks$Consensus)
ConsensusPeaks <- dba.peakset(ATACConsensus, bRetrieve=TRUE)
ConsensusPeaks
save(ConsensusPeaks,file="PSEN1.PSEN2.APP.ATACseq.ConsensusPeaks.Rdata")
# Annotate
anno.ConsensusPeaks <- annotatePeak(ConsensusPeaks, 
                                    tssRegion=c(-1000,500), 
                                    TxDb=txdb.hg38ens, 
                                    level="gene", 
                                    annoDb="org.Hs.eg.db",
                                    overlap="TSS")
anno.ConsensusPeaks.ranges <- as_granges(anno.ConsensusPeaks)
# Define promoter and enhancer regions
ConsensusPeaks.promoter <- anno.ConsensusPeaks.ranges[anno.ConsensusPeaks.ranges$annotation == "Promoter"]
ConsensusPeaks.promoter
ConsensusPeaks.nonpromoter <- anno.ConsensusPeaks.ranges[anno.ConsensusPeaks.ranges$annotation != "Promoter"]
ConsensusPeaks.nonpromoter
ConsensusPeaks.enhancer.inner <- ConsensusPeaks.nonpromoter %>%
  join_overlap_inner(PEREGRINE.brain.ranges)
# center scale function
center_rowmeans <- function(x) {
  xcenter = rowMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}
center_scale <- function(x) {
  scale(x, scale = FALSE)
}

# Extract the ATAC-seq voom counts
ATAC.log.counts <- v$E
ATAC.peak.regions <- ATACcounts.df[,c(1:3)]
ATAC.peak.regions$start.end <- paste(ATAC.peak.regions$start, ATAC.peak.regions$end, sep="-")
ATAC.peak.regions$combined <- paste(ATAC.peak.regions$seqnames, ATAC.peak.regions$start.end, sep=":")
ATAC.peak.regions.combined <- as.data.frame(ATAC.peak.regions$combined)
colnames(ATAC.peak.regions.combined) <- "loc"
ATAC.log.counts.df <- as.data.frame(ATAC.log.counts)
ATAC.log.counts.df$NDC <- rowMeans(ATAC.log.counts[,1:3])
ATAC.log.counts.df$V717I <- rowMeans(ATAC.log.counts[,4:6])
ATAC.log.counts.df$A79V <- rowMeans(ATAC.log.counts[,7:9])
ATAC.log.counts.df$N141 <- rowMeans(ATAC.log.counts[,10:12])
ATAC.log.counts.avg <- ATAC.log.counts.df[,c(13:16)]
ATAC.log.counts.meancentered <- center_scale(ATAC.log.counts.avg)
ATAC.maelstrom <- dplyr::bind_cols(ATAC.peak.regions.combined,as.data.frame(ATAC.log.counts.meancentered))
head(ATAC.maelstrom)
# Use only the consensus peak regions
ATAC.summits.coverage <- dplyr::bind_cols(ATACcounts.df[,c(1:3)],as.data.frame(ATAC.log.counts.meancentered))
ATAC.summits.coverage.ranges <- as_granges(ATAC.summits.coverage)
ConsensusPeaks.df <- as.data.frame(ConsensusPeaks)
ConsensusPeaks.3col.df <- ConsensusPeaks.df[,c(1:3)]
ConsensusPeaks.3col <- as_granges(ConsensusPeaks.3col.df)
ATAC.consensus.coverage.ranges <-  ATAC.summits.coverage.ranges %>%
  join_overlap_inner(ConsensusPeaks.3col)
ATAC.consensus.peaks <- as.data.frame(ATAC.consensus.coverage.ranges)[,c(1:3)]
ATAC.consensus.peaks$start.end <- paste(ATAC.consensus.peaks$start, ATAC.consensus.peaks$end, sep="-")
ATAC.consensus.peaks$combined <- paste(ATAC.consensus.peaks$seqnames, ATAC.consensus.peaks$start.end, sep=":")
ATAC.consensus.peaks.combined <- as.data.frame(ATAC.consensus.peaks$combined)
colnames(ATAC.consensus.peaks.combined) <- "loc"
ATAC.consensus.maelstrom <- dplyr::bind_cols(ATAC.consensus.peaks.combined,as.data.frame(ATAC.consensus.coverage.ranges)[,c(6:9)])
ATAC.consensus.maelstrom
# Export data matrix for GimmeMotifs
write.table(ATAC.consensus.maelstrom, file = "PSEN1.PSEN2.APP.consensusATAC.GimmeMotifs.Maelstrom.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# Subset the promoterpeaks
ConsensusPeaks.promoter.df <- as.data.frame(ConsensusPeaks.promoter)
ConsensusPeaks.promoter.3col.df <- ConsensusPeaks.promoter.df[,c(1:3)]
ConsensusPeaks.promoter.3col <- as_granges(ConsensusPeaks.promoter.3col.df)
ATAC.consensus.promoter.coverage.ranges <-  ATAC.summits.coverage.ranges %>%
  join_overlap_inner(ConsensusPeaks.promoter.3col)
ATAC.peaks.promoter <- as.data.frame(ATAC.consensus.promoter.coverage.ranges)[,c(1:3)]
ATAC.peaks.promoter$start.end <- paste(ATAC.peaks.promoter$start, ATAC.peaks.promoter$end, sep="-")
ATAC.peaks.promoter$combined <- paste(ATAC.peaks.promoter$seqnames, ATAC.peaks.promoter$start.end, sep=":")
ATAC.peaks.promoter.combined <- as.data.frame(ATAC.peaks.promoter$combined)
colnames(ATAC.peaks.promoter.combined) <- "loc"
ATAC.consensus.promoter.maelstrom <- dplyr::bind_cols(ATAC.peaks.promoter.combined,as.data.frame(ATAC.consensus.promoter.coverage.ranges)[,c(6:9)])
ATAC.consensus.promoter.maelstrom
# Export promoter data matrix for GimmeMotifs
write.table(ATAC.consensus.promoter.maelstrom, file = "PSEN1.PSEN2.APP.consensus.promoter.Maelstrom.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# enhancer.inner 
ConsensusPeaks.enhancer.inner.df <- as.data.frame(unique(ConsensusPeaks.enhancer.inner))
ConsensusPeaks.enhancer.inner.3col.df <- ConsensusPeaks.enhancer.inner.df[,c(1:3)]
ConsensusPeaks.enhancer.inner.3col <- as_granges(ConsensusPeaks.enhancer.inner.3col.df)
ATAC.consensus.enhancer.inner.coverage.ranges <-  ATAC.summits.coverage.ranges %>%
  join_overlap_inner(ConsensusPeaks.enhancer.inner.3col)
ATAC.peaks.enhancer.inner <- as.data.frame(ATAC.consensus.enhancer.inner.coverage.ranges)[,c(1:3)]
ATAC.peaks.enhancer.inner$start.end <- paste(ATAC.peaks.enhancer.inner$start, ATAC.peaks.enhancer.inner$end, sep="-")
ATAC.peaks.enhancer.inner$combined <- paste(ATAC.peaks.enhancer.inner$seqnames, ATAC.peaks.enhancer.inner$start.end, sep=":")
ATAC.peaks.enhancer.inner.combined <- as.data.frame(ATAC.peaks.enhancer.inner$combined)
colnames(ATAC.peaks.enhancer.inner.combined) <- "loc"
ATAC.consensus.enhancer.inner.maelstrom <- dplyr::bind_cols(ATAC.peaks.enhancer.inner.combined,as.data.frame(ATAC.consensus.enhancer.inner.coverage.ranges)[,c(6:9)])
ATAC.consensus.enhancer.inner.maelstrom
# Export enhancer data matrix for GimmeMotifs
write.table(ATAC.consensus.enhancer.inner.maelstrom, file = "PSEN1.PSEN2.APP.consensus.enhancer.inner.Maelstrom.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# save results
write.table(ATAC.maelstrom, file = "PSEN1.PSEN2.APP.GimmeMotifs.Maelstrom.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
rtracklayer::export.bed((ATAC.peak.regions[,c(1:3)]), con = "PSEN1.PSEN2.APP.allATACpeaks.bed", format = "bed")
