#RNA-seq/ATAC-seq data integration
#Valdez et. al 2023 Molecular Psychiatry Submission
#load packages
library("DiffBind")
library("ChIPseeker")
library("rtracklayer")
library("plyranges")
library("edgeR")
library("limma")
library("msigdbr")
library("tmod")

#load hg38
txdb.hg38ens <- loadDb("~/gmtfiles/txdb.hg38ens.sqlite")
#load the RNA-seq data
load("~/RNA-seq/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")
#load the ATAC-seq data
load("~/ATAC-seq/ATACSeq.report.V717I.Rdata")
load("~/ATAC-seq/ATACSeq.report.A79V.Rdata")
load("~/ATAC-seq/ATACSeq.report.N141I.Rdata")

#
anno.V717I <- annotatePeak(ATACSeq.report.V717I, 
                           tssRegion=c(-1500,500), 
                           TxDb=txdb.hg38ens, 
                           level="gene", 
                           annoDb="org.Hs.eg.db",
                           overlap="TSS")
anno.V717I.ranges <- as_granges(anno.V717I)
# V717I
anno.V717I.ranges <- anno.V717I.ranges %>%
  mutate(
#    da_log2FC = Fold,
 #   da_padj = FDR,
    ENSEMBL = geneId 
  )
# Select key columns and convert back to granges
dATACpeaks.V717I <- as_granges(as.data.frame(anno.V717I.ranges)[,c(1:5,12,9,11,23,21)])
# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.V717I) <- dATACpeaks.V717I$ENSEMBL
dATACpeaks.V717I
# View the significant differential ATAC peaks
dATACpeaks.V717I.ranges <- dATACpeaks.V717I
# Use the treat method for differential gene expression
sigGenes.V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)
sigGenes.V717I.filt <- sigGenes.V717I %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(abs(logFC) > 0)
length(unique(sigGenes.V717I.filt$ENSEMBL))
# Isolate the key attributes and rename the columns
sigGenes.V717I <- sigGenes.V717I %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.V717I)
# Remove duplicated ENSEMBL terms
sigGenes.V717I = subset(sigGenes.V717I, ENSEMBL != "" )
sigGenes.V717I = subset(sigGenes.V717I, ! duplicated(ENSEMBL))
dim(sigGenes.V717I)
head(sigGenes.V717I)
# Match the genes from gene expresssion data to the differentially accessible ATAC peaks
dATACpeaks.V717I.w <- dATACpeaks.V717I.ranges[na.omit(match(rownames(sigGenes.V717I), names(dATACpeaks.V717I.ranges)))]
# Add the differential expression statistics of genes ID'd by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.V717I.w) <- sigGenes.V717I[match(names(dATACpeaks.V717I.w), rownames(sigGenes.V717I)),]
# View the intergrated results
dATACpeaks.V717I.w
# Overlap the dATAC with dE peak regions
dATACpeaks.V717I.final <- dATACpeaks.V717I %>%
  join_overlap_left(dATACpeaks.V717I.w)
# Filter for both differential gene expression and differential accessibility
dATACpeaks.V717I.DEG  <- dATACpeaks.V717I.final  %>%
  filter(de_padj < 0.05) %>%
  filter(abs(de_log2FC) > 0) %>%
  filter(da_padj < 0.05) %>%
  filter(abs(da_log2FC) > 0)
# View the final granges object of dATAC peaks with differential gene expression
dATACpeaks.V717I.DEG 
# Isolate the key statistics of each
dATACpeaks.V717I.DEG <- dATACpeaks.V717I.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.V717I.DEG
save(dATACpeaks.V717I.DEG,file="~/RNA_ATAC/dATACpeaks.V717I.DEG.Rdata")
write.table(dATACpeaks.V717I.DEG,file="~/RNA_ATAC/dATACpeaks.V717I.DEG.tsv",sep="\t")
#
anno.V717I.ranges <- anno.V717I.ranges %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )
# Select key columns and convert back to granges
dATACpeaks.V717I <- as_granges(as.data.frame(anno.V717I.ranges)[,c(1:8,12:25)])
# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.V717I) <- dATACpeaks.V717I$ENSEMBL
dATACpeaks.V717I
# Isolate the peaks with only ENSEMBL IDs as metadata for matching with gene expression data
dATACpeaks.V717I.ranges <- as_granges(anno.V717I[,c(1:5,25)])
# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.V717I.ranges) <- dATACpeaks.V717I$ENSEMBL
# View the significant differential ATAC peaks
dATACpeaks.V717I.ranges <- dATACpeaks.V717I
# Use the treat method for differential gene expression
sigGenes.V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

# Isolate the key attributes and rename the columns
sigGenes.V717I <- sigGenes.V717I %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.V717I)
# Remove duplicated ENSEMBL terms
sigGenes.V717I = subset(sigGenes.V717I, ENSEMBL != "" )
sigGenes.V717I = subset(sigGenes.V717I, ! duplicated(ENSEMBL))
dim(sigGenes.V717I)
head(sigGenes.V717I)
# Match the genes from gene expresssion data to the differentially accessible ATAC peaks
dATACpeaks.V717I.w <- dATACpeaks.V717I.ranges[na.omit(match(rownames(sigGenes.V717I), names(dATACpeaks.V717I.ranges)))]
# Add the differential expression statistics of genes ID'd by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.V717I.w) <- sigGenes.V717I[match(names(dATACpeaks.V717I.w), rownames(sigGenes.V717I)),]
# View the intergrated results
dATACpeaks.V717I.w
# Overlap the dATAC with dE peak regions
dATACpeaks.V717I.final <- dATACpeaks.V717I %>%
  join_overlap_left(dATACpeaks.V717I.w)
# Filter for both differential gene expression and differential accessibility
dATACpeaks.V717I.DEG  <- dATACpeaks.V717I.final  %>%
  filter(de_padj < 0.05) %>%
  filter(abs(de_log2FC) > 0) %>%
  filter(da_padj < 0.05) %>%
  filter(abs(da_log2FC) > 0)
# View the final granges object of dATAC peaks with differential gene expression
dATACpeaks.V717I.DEG 
# Isolate the key statistics of each
dATACpeaks.V717I.DEG <- dATACpeaks.V717I.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.V717I.DEG
save(dATACpeaks.V717I.DEG,file="~/RNA_ATAC/dATACpeaks.V717I.voomDEG.Rdata")
write.table(dATACpeaks.V717I.DEG,file="~/RNA_ATAC/dATACpeaks.V717I.voomDEG.tsv",sep="\t")
# A79V
anno.A79V.ranges <- anno.A79V.ranges %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )
# Select key columns and convert back to granges
dATACpeaks.A79V <- as_granges(as.data.frame(anno.A79V.ranges)[,c(1:8,12:25)])
# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.A79V) <- dATACpeaks.A79V$ENSEMBL
dATACpeaks.A79V
# Isolate the peaks with only ENSEMBL IDs as metadata for matching with gene expression data
dATACpeaks.A79V.ranges <- as_granges(anno.A79V[,c(1:5,25)])
# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.A79V.ranges) <- dATACpeaks.A79V$ENSEMBL
# View the significant differential ATAC peaks
dATACpeaks.A79V.ranges <- dATACpeaks.A79V
# Use the treat method for differential gene expression
sigGenes.A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)
sigGenes.A79V.filt <- sigGenes.A79V %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(abs(logFC) > 0)
length(unique(sigGenes.A79V.filt$ENSEMBL))
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
# Add the differential expression statistics of genes ID'd by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.A79V.w) <- sigGenes.A79V[match(names(dATACpeaks.A79V.w), rownames(sigGenes.A79V)),]
# View the intergrated results
dATACpeaks.A79V.w
# Overlap the dATAC with dE peak regions
dATACpeaks.A79V.final <- dATACpeaks.A79V %>%
  join_overlap_left(dATACpeaks.A79V.w)
# Filter for both differential gene expression and differential accessibility
dATACpeaks.A79V.DEG  <- dATACpeaks.A79V.final  %>%
  filter(de_padj < 0.05) %>%
  filter(abs(de_log2FC) > 0) %>%
  filter(da_padj < 0.05) %>%
  filter(abs(da_log2FC) > 0)
# View the final granges object of dATAC peaks with differential gene expression
dATACpeaks.A79V.DEG 
# Isolate the key statistics of each
dATACpeaks.A79V.DEG <- dATACpeaks.A79V.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.A79V.DEG
save(dATACpeaks.A79V.DEG,file="~/RNA_ATAC/dATACpeaks.A79V.voomDEG.Rdata")
write.table(dATACpeaks.A79V.DEG,file="~/RNA_ATAC/dATACpeaks.A79V.voomDEG.tsv",sep="\t")
#
anno.N141I.ranges <- anno.N141I.ranges %>%
  mutate(
    da_log2FC = Fold,
    da_padj = FDR,
    ENSEMBL = geneId 
  )
# Select key columns and convert back to granges
dATACpeaks.N141I <- as_granges(as.data.frame(anno.N141I.ranges)[,c(1:8,12:25)])
# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.N141I) <- dATACpeaks.N141I$ENSEMBL
dATACpeaks.N141I
# Isolate the peaks with only ENSEMBL IDs as metadata for matching with gene expression data
dATACpeaks.N141I.ranges <- as_granges(anno.N141I[,c(1:5,25)])
# Add ENSEMBL IDs as row names of the granges object
names(dATACpeaks.N141I.ranges) <- dATACpeaks.N141I$ENSEMBL
# View the significant differential ATAC peaks
dATACpeaks.N141I.ranges <- dATACpeaks.N141I
# Use the treat method for differential gene expression
sigGenes.N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)
sigGenes.N141I.filt <- sigGenes.N141I %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::filter(abs(logFC) > 0)
length(unique(sigGenes.N141I.filt$ENSEMBL))
# Isolate the key attributes and rename the columns
sigGenes.N141I <- sigGenes.N141I %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.N141I)
# Remove duplicated ENSEMBL terms
sigGenes.N141I = subset(sigGenes.N141I, ENSEMBL != "" )
sigGenes.N141I = subset(sigGenes.N141I, ! duplicated(ENSEMBL))
dim(sigGenes.N141I)
head(sigGenes.N141I)
# Match the genes from gene expresssion data to the differentially accessible ATAC peaks
dATACpeaks.N141I.w <- dATACpeaks.N141I.ranges[na.omit(match(rownames(sigGenes.N141I), names(dATACpeaks.N141I.ranges)))]
# Add the differential expression statistics of genes ID'd by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.N141I.w) <- sigGenes.N141I[match(names(dATACpeaks.N141I.w), rownames(sigGenes.N141I)),]
# View the intergrated results
dATACpeaks.N141I.w
# Overlap the dATAC with dE peak regions
dATACpeaks.N141I.final <- dATACpeaks.N141I %>%
  join_overlap_left(dATACpeaks.N141I.w)
# Filter for both differential gene expression and differential accessibility
dATACpeaks.N141I.DEG  <- dATACpeaks.N141I.final  %>%
  filter(de_padj < 0.05) %>%
  filter(abs(de_log2FC) > 0) %>%
  filter(da_padj < 0.05) %>%
  filter(abs(da_log2FC) > 0)
# View the final granges object of dATAC peaks with differential gene expression
dATACpeaks.N141I.DEG 
# Isolate the key statistics of each
dATACpeaks.N141I.DEG <- dATACpeaks.N141I.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.N141I.DEG
save(dATACpeaks.N141I.DEG,file="~/RNA_ATAC/dATACpeaks.N141I.voomDEG.Rdata")
write.table(dATACpeaks.N141I.DEG,file="~/RNA_ATAC/dATACpeaks.N141I.voomDEG.tsv",sep="\t")
#
# Create VennEuler overlap diagrams
library("venneuler")
library("VennDiagram")
library("eulerr")
# V717I
venn.plot.V717I <- draw.pairwise.venn(area1 = length(unique(sigGenes.V717I.filt$ENSEMBL)),
                                area2 = length(unique(names(dATACpeaks.V717I.ranges))),
                                cross.area = length(unique(names(dATACpeaks.V717I.DEG))),
                                euler.d = TRUE,
                                scaled = TRUE)
length(unique(sigGenes.V717I.filt$ENSEMBL))
length(unique(names(dATACpeaks.V717I.ranges)))
length(unique(names(dATACpeaks.V717I.DEG)))
venn.euler.V717I <- eulerr::euler(c(A = length(unique(sigGenes.V717I.filt$ENSEMBL)),
                                      B = length(unique(names(dATACpeaks.V717I.ranges))),
                                      "A&B" = length(unique(names(dATACpeaks.V717I.DEG)))))
plot(venn.euler.V717I)
# A79V
dev.off()
venn.plot.A79V <- draw.pairwise.venn(area1 = length(unique(sigGenes.A79V.filt$ENSEMBL)),
                                      area2 = length(unique(names(dATACpeaks.A79V.ranges))),
                                      cross.area = length(unique(names(dATACpeaks.A79V.DEG))),
                                      euler.d = TRUE,
                                      scaled = TRUE)
length(unique(sigGenes.A79V.filt$ENSEMBL))
length(unique(names(dATACpeaks.A79V.ranges)))
length(unique(names(dATACpeaks.A79V.DEG)))
venn.euler.A79V <- eulerr::euler(c(A = length(unique(sigGenes.A79V.filt$ENSEMBL)),
                                    B = length(unique(names(dATACpeaks.A79V.ranges))),
                                    "A&B" = length(unique(names(dATACpeaks.A79V.DEG)))))
plot(venn.euler.A79V)
# N141I
dev.off()
venn.plot.N141I <- draw.pairwise.venn(area1 = length(unique(sigGenes.N141I.filt$ENSEMBL)),
                                      area2 = length(unique(names(dATACpeaks.N141I.ranges))),
                                      cross.area = length(unique(names(dATACpeaks.N141I.DEG))),
                                      euler.d = TRUE,
                                      scaled = TRUE)
length(unique(sigGenes.N141I.filt$ENSEMBL))
length(unique(names(dATACpeaks.N141I.ranges)))
length(unique(names(dATACpeaks.N141I.DEG)))
venn.euler.N141I <- eulerr::euler(c(A = length(unique(sigGenes.N141I.filt$ENSEMBL)),
                                   B = length(unique(names(dATACpeaks.N141I.ranges))),
                                   "A&B" = length(unique(names(dATACpeaks.N141I.DEG)))))
plot(venn.euler.N141I)
#
#
dATACpeaks.A79V.DEG.down  <- dATACpeaks.A79V.DEG  %>%
  filter(de_log2FC < 0) %>%
  filter(da_log2FC < 0)
dATACpeaks.A79V.DEG.up  <- dATACpeaks.A79V.DEG  %>%
  filter(de_log2FC > 0) %>%
  filter(da_log2FC > 0)

#####################################################################
# Run hypergeometric enrichment on integrated differential RNA-seq/ATAC-seq genes
load("~/ATAC-seq/PSEN1.PSEN2.APP.ATACSeq.report.A79V.full.Rdata")
# Get the background ATAC-seq genes
A79V.all.granges <- as_granges(ATACSeq.report.A79V_full)
anno.A79V.all.granges <- annotatePeak(A79V.all.granges, 
                                      tssRegion=c(-1000,500), 
                                      TxDb=txdb.hg38ens, 
                                      level="gene", 
                                      annoDb="org.Hs.eg.db",
                                      overlap="TSS")
anno.A79V.all.granges <- as_granges(anno.A79V.all.granges)
dATAC.genes.background <- unique(as.character(anno.A79V.all.granges$SYMBOL))
#
#set V717I dATAC/dEG genes
dATACpeaks.V717I.DEG.down  <- dATACpeaks.V717I.DEG  %>%
  filter(de_log2FC < 0) %>%
  filter(da_log2FC < 0)
dATACpeaks.V717I.DEG.up  <- dATACpeaks.V717I.DEG  %>%
  filter(de_log2FC > 0) %>%
  filter(da_log2FC > 0)
# set genes
down_genes <- unique(as.character(dATACpeaks.V717I.DEG.down$SYMBOL))
up_genes <- unique(as.character(dATACpeaks.V717I.DEG.up$SYMBOL))
# set background genes
all_genes <- union(unique(as.character(y$genes$SYMBOL)),dATAC.genes.background)
# load the TF databasets
load("~/gmtfiles/msetECC.Rdata")
load("~/gmtfiles/msetReMap.Rdata")
#Import the latest MSigDB XML File -> download 7.4 version
msig <- tmodImportMSigDB(file.path("~/msigdb", "msigdb_v7.4.xml"))
#Select the GO BP gene sets (you can select whichever geneset db you would like)
Hallmark.sel <- msig$MODULES$Category == "H"
GOBP.exp.sel <- msig$MODULES$Subcategory == "GO:BP"
# ReMap TF database
tmodHG_ReMap_down <- tmodHGtest(down_genes, all_genes, mset = msetReMap)
head(tmodHG_ReMap_down, 20)
tmodHG_ReMap_up <- tmodHGtest(up_genes, all_genes, mset = msetReMap)
head(tmodHG_ReMap_up, 20)
#ECC TF database
tmodHG_ECC_down <- tmodHGtest(down_genes, all_genes, mset = msetECC)
head(tmodHG_ECC_down, 20)
tmodHG_ECC_up <- tmodHGtest(up_genes, all_genes, mset = msetECC)
head(tmodHG_ECC_up, 20)
#GOBP
tmodHG_GOBP.exp_down <- tmodHGtest(down_genes, all_genes, mset = msig[GOBP.exp.sel])
head(tmodHG_GOBP.exp_down, 20)
tmodHG_GOBP.exp_up <- tmodHGtest(up_genes, all_genes, mset = msig[GOBP.exp.sel])
head(tmodHG_GOBP.exp_up, 20)
#Hallmark
tmodHG_Hallmark_down <- tmodHGtest(down_genes, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark_down, 20)
tmodHG_Hallmark_up <- tmodHGtest(up_genes, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark_up, 20)
write.table(tmodHG_ReMap_down,file="~/RNA_ATAC/DEG.dATAC.V717I.ReMap.down.tsv",sep="\t")
write.table(tmodHG_ReMap_up,file="~/RNA_ATAC/DEG.dATAC.V717I.ReMap.up.tsv",sep="\t")
write.table(tmodHG_ECC_down,file="~/RNA_ATAC/DEG.dATAC.V717I.ECC.down.tsv",sep="\t")
write.table(tmodHG_ECC_up,file="~/RNA_ATAC/DEG.dATAC.V717I.ECC.up.tsv",sep="\t")
write.table(tmodHG_GOBP.exp_down,file="~/RNA_ATAC/DEG.dATAC.V717I.GOBPexp.down.tsv",sep="\t")
write.table(tmodHG_GOBP.exp_up,file="~/RNA_ATAC/DEG.dATAC.V717I.GOBPexp.up.tsv",sep="\t")
write.table(tmodHG_Hallmark_down,file="~/RNA_ATAC/DEG.dATAC.V717I.Hallmark.down.tsv",sep="\t")
write.table(tmodHG_Hallmark_up,file="~/RNA_ATAC/DEG.dATAC.V717I.Hallmark.up.tsv",sep="\t")
#####################################################################
# set A79V dATAC/dEG genes
down_genes <- unique(as.character(dATACpeaks.A79V.DEG.down$SYMBOL))
up_genes <- unique(as.character(dATACpeaks.A79V.DEG.up$SYMBOL))
#GOBP.exp
tmodHG_ReMap_down <- tmodHGtest(down_genes, all_genes, mset = msetReMap)
head(tmodHG_ReMap_down, 20)
tmodHG_ReMap_up <- tmodHGtest(up_genes, all_genes, mset = msetReMap)
head(tmodHG_ReMap_up, 20)
#ECC
tmodHG_ECC_down <- tmodHGtest(down_genes, all_genes, mset = msetECC)
head(tmodHG_ECC_down, 20)
tmodHG_ECC_up <- tmodHGtest(up_genes, all_genes, mset = msetECC)
head(tmodHG_ECC_up, 20)
#GOBP
tmodHG_GOBP.exp_down <- tmodHGtest(down_genes, all_genes, mset = msig[GOBP.exp.sel])
head(tmodHG_GOBP.exp_down, 20)
tmodHG_GOBP.exp_up <- tmodHGtest(up_genes, all_genes, mset = msig[GOBP.exp.sel])
head(tmodHG_GOBP.exp_up, 20)
#Hallmark
tmodHG_Hallmark_down <- tmodHGtest(down_genes, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark_down, 20)
tmodHG_Hallmark_up <- tmodHGtest(up_genes, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark_up, 20)
write.table(tmodHG_ReMap_down,file="~/RNA_ATAC/DEG.dATAC.A79V.ReMap.down.tsv",sep="\t")
write.table(tmodHG_ReMap_up,file="~/RNA_ATAC/DEG.dATAC.A79V.ReMap.up.tsv",sep="\t")
write.table(tmodHG_ECC_down,file="~/RNA_ATAC/DEG.dATAC.A79V.ECC.down.tsv",sep="\t")
write.table(tmodHG_ECC_up,file="~/RNA_ATAC/DEG.dATAC.A79V.ECC.up.tsv",sep="\t")
write.table(tmodHG_GOBP.exp_down,file="~/RNA_ATAC/DEG.dATAC.A79V.GOBPexp.down.tsv",sep="\t")
write.table(tmodHG_GOBP.exp_up,file="~/RNA_ATAC/DEG.dATAC.A79V.GOBPexp.up.tsv",sep="\t")
write.table(tmodHG_ReMap_down,file="~/RNA_ATAC/DEG.dATAC.A79V.ReMap.down.tsv",sep="\t")
write.table(tmodHG_ReMap_up,file="~/RNA_ATAC/DEG.dATAC.A79V.ReMap.up.tsv",sep="\t")
#####################################################################
#set N141I dATAC/dEG
dATACpeaks.N141I.DEG.down  <- dATACpeaks.N141I.DEG  %>%
  filter(de_log2FC < 0) %>%
  filter(da_log2FC < 0)
dATACpeaks.N141I.DEG.up  <- dATACpeaks.N141I.DEG  %>%
  filter(de_log2FC > 0) %>%
  filter(da_log2FC > 0)
down_genes <- unique(as.character(dATACpeaks.N141I.DEG.down$SYMBOL))
up_genes <- unique(as.character(dATACpeaks.N141I.DEG.up$SYMBOL))
#GOBP.exp
tmodHG_ReMap_down <- tmodHGtest(down_genes, all_genes, mset = msetReMap)
head(tmodHG_ReMap_down, 20)
tmodHG_ReMap_up <- tmodHGtest(up_genes, all_genes, mset = msetReMap)
head(tmodHG_ReMap_up, 20)
#ECC
tmodHG_ECC_down <- tmodHGtest(down_genes, all_genes, mset = msetECC)
head(tmodHG_ECC_down, 20)
tmodHG_ECC_up <- tmodHGtest(up_genes, all_genes, mset = msetECC)
head(tmodHG_ECC_up, 20)
#GOBP
tmodHG_GOBP.exp_down <- tmodHGtest(down_genes, all_genes, mset = msig[GOBP.exp.sel])
head(tmodHG_GOBP.exp_down, 20)
tmodHG_GOBP.exp_up <- tmodHGtest(up_genes, all_genes, mset = msig[GOBP.exp.sel])
head(tmodHG_GOBP.exp_up, 20)
#Hallmark
tmodHG_Hallmark_down <- tmodHGtest(down_genes, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark_down, 20)
tmodHG_Hallmark_up <- tmodHGtest(up_genes, all_genes, mset = msig[Hallmark.sel])
head(tmodHG_Hallmark_up, 20)
write.table(tmodHG_ReMap_down,file="~/RNA_ATAC/DEG.dATAC.N141I.ReMap.down.tsv",sep="\t")
write.table(tmodHG_ReMap_up,file="~/RNA_ATAC/DEG.dATAC.N141I.ReMap.up.tsv",sep="\t")
write.table(tmodHG_ECC_down,file="~/RNA_ATAC/DEG.dATAC.N141I.ECC.down.tsv",sep="\t")
write.table(tmodHG_ECC_up,file="~/RNA_ATAC/DEG.dATAC.N141I.ECC.up.tsv",sep="\t")
write.table(tmodHG_GOBP.exp_down,file="~/RNA_ATAC/DEG.dATAC.N141I.GOBPexp.down.tsv",sep="\t")
write.table(tmodHG_GOBP.exp_up,file="~/RNA_ATAC/DEG.dATAC.N141I.GOBPexp.up.tsv",sep="\t")
write.table(tmodHG_ReMap_down,file="~/RNA_ATAC/DEG.dATAC.N141I.ReMap.down.tsv",sep="\t")
write.table(tmodHG_ReMap_up,file="~/RNA_ATAC/DEG.dATAC.N141I.ReMap.up.tsv",sep="\t")
