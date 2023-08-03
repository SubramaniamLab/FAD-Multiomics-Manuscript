#ATAC-seq intepareto analysis 
#Valdes et. al 2023 Molecular Psychiatry Submission
#load packages
library("intePareto")
library("biomaRt")
library("tmod")
library("msigdbr")
# set working dirECCtory
base_dir <- getwd()
# Run intePareto for each comparison: NDC-A79V (PSEN1), NDC-N141I (PSEN2), and NDC-V717I (APP)
# Example is for NDC-A79V
# 3 column txt file of the RNA abundance.tsv kallisto output of the two comparison group (colnames: SRR, condition, files)
rnaMeta <- read.table(file.path(base_dir, "NDC_A79V_intePareto_RNA.txt"), header = TRUE, stringsAsFactors=FALSE)
# 3 column txt file of the promoter-only bam files of the two comparison group (colnames: SRR, condition, files)
# subset bams for GRCh38 promoters using bedtools in command line (repeat for all replicates)
# First need to do bedtools slop to expand TSS to putative promoter
#Use FANTOM5 promoters
#awk -F'\t' 'OFS="\t" {print $1,$7,$8}' FANTOM5.liftover.CAGE.promoter.hg38.bed > FANTOM5.liftover.CAGE.TSS.hg38.bed
#bedtools slop \
#-i FANTOM5.liftover.CAGE.TSS.hg38.bed \
#-l 1000 \
#-r 500 \
#-g /home/andrewbcaldwell/rgtdata/hg38/chrom.sizes.hg38 \
#> FANTOM5.liftover.CAGE.TSS.promoter.hg38.intePareto.bed
#bedtools intersECCt \
#-abam 27_4-NDC1-1_S1_combined.bb.mapped.sorted.bam \
#-b FANTOM5.liftover.CAGE.TSS.promoter.hg38.intePareto.bed > 27_4-NDC1-1_S1.promoter.bam
#done
# Make 3 column txt file of the promoter bam files of the two comparison group (colnames: SRR, condition, files)
patacMeta <- read.table(file.path(base_dir, "NDC_A79V_intePareto_promoterATAC.txt"), header = TRUE, stringsAsFactors=FALSE)
# subset bams for PEREGRINE ehancers using bedtools in command line (repeat for all replicates)
#bedtools intersECCt \
#-abam 27_4-NDC1-1_S1_combined.combined.bb.mapped.sorted.bam \
#-b ~/gmtfiles/PEREGRINE.brain.enhancer.bed > 27_4-NDC1-1_S1_combined.PEREGRINE.enhancer.bam
#done
# Make 3 column txt file of the enhancer bam files of the two comparison group (colnames: SRR, condition, files)
eatacMeta <- read.table(file.path(base_dir, "NDC_A79V_intePareto_enhancerATAC.txt"), header = TRUE, stringsAsFactors=FALSE)
# Run intePareto
#promoter
res.p <- doMatch(rnaMeta = rnaMeta,
               chipMeta = patacMeta,
               region = "promoter",
               method = "highest",
               ensemblDataset = "hsapiens_gene_ensembl",
               host = "jan2020.archive.ensembl.org",
               fragLength = 150, 
               promoter.length = 3000)
#enhancer
res.e <- doMatch(rnaMeta = rnaMeta,
                 chipMeta = eatacMeta,
                 region = "genebody",
                 method = "highest",
                 ensemblDataset = "hsapiens_gene_ensembl",
                 host = "jan2020.archive.ensembl.org")
# Merge
res.p$matched.data <- merge(res.p$matched.data,
                            res.e$matched.data)
res.p$res.chip <- merge(res.p$res.chip,
                        res.e$res.chip)
res <- res.p
# integration
res.integration <- doIntegration(res,
                                 ref = "NDC",
                                 type = "apeglm",
                                 apeAdapt = FALSE)
# define fronts
nr.fronts <- length(res.integration$RNAseq.log2FoldChange)
# Define objECCtives
objECCtive <-data.frame(mark =c("z.pATAC","z.eATAC"),
                       obj=c("max","max"),
                       stringsAsFactors=FALSE)
# Pareto optimization
res.final <-doPareto(df_final = res.integration,
                     objECCtive = objECCtive,
                     nr.fronts = 50000)
# Final results
res.final
# save results
write.csv(res.final, file = "intePareto.A79V.results.csv") 
# Filter results
res.final$SYMBOL <- row.names(res.final)
res.final.df <- as.data.frame(res.final)
res.final.dt = subset(res.final.df, SYMBOL != "" )
res.final.dt = subset(res.final.df, ! duplicated(SYMBOL))
# Tidy up and order by t metric
load("~/RNA-seq/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")
limma_A79V <- topTable(efit, coef="A79VvNDC", adjust.method = "BH", sort.by ="t", n = Inf) 
limma_A79V$SYMBOL <- limma_A79V$SYMBOL
ranks_intePareto <- res.final.dt
row.names.remove <- setdiff(ranks_intePareto$SYMBOL,limma_A79V$SYMBOL)
ranks_intePareto <- ranks_intePareto[!(ranks_intePareto$SYMBOL %in% row.names.remove), ]
ranks_intePareto <- subset(ranks_intePareto, ! duplicated(ranks_intePareto$SYMBOL))
ranks <- deframe(ranks_intePareto)
# Run CERNO enrichment on the ordered intePareto results
#Import the latest MSigDB XML File -> download 7.4 version
msig <- tmodImportMSigDB(file.path("~/msigdb", "msigdb_v7.4.xml"))
#SelECCt the GO BP gene sets (you can selECCt whichever geneset db you would like)
Hallmark.sel <- msig$MODULES$Category == "H"
GOBP.exp.sel <- msig$MODULES$Subcategory == "GO:BP"
# Hallmark
CERNO.Hallmark <- tmodCERNOtest(
  ranks_intePareto$SYMBOL,
  modules = NULL,
  qval = 0.05,
  order.by = "pval",
  filter = FALSE,
  mset = msig[Hallmark.sel],
  cols = "Title")
CERNO.Hallmark
write.csv(CERNO.Hallmark, file = "intePareto.CERNO.A79V.Hallmark.csv") 
# GOBP
CERNO.GOBP <- tmodCERNOtest(
  ranks_intePareto$SYMBOL,
  modules = NULL,
  qval = 0.05,
  order.by = "pval",
  filter = FALSE,
  mset = msig[GOBP.sel],
  cols = "Title")
CERNO.GOBP
write.csv(CERNO.GOBP, file = "intePareto.CERNO.A79V.GOBP.csv") 
# load the TF databasets
load("~/gmtfiles/msetECCC.Rdata")
load("~/gmtfiles/msetReMap.Rdata")
# ReMap
CERNO.ReMap <- tmodCERNOtest(
  ranks_intePareto$SYMBOL,
  modules = NULL,
  qval = 0.05,
  order.by = "pval",
  filter = FALSE,
  mset = msetReMap,
  cols = "Title")
CERNO.ReMap
write.csv(CERNO.ReMap, file = "intePareto.CERNO.A79V.ReMap.csv") 
# ECC
CERNO.ECC <- tmodCERNOtest(
  ranks_intePareto$SYMBOL,
  modules = NULL,
  qval = 0.05,
  order.by = "pval",
  filter = FALSE,
  mset = msetECC,
  cols = "Title")
CERNO.ECC
write.csv(CERNO.ECC, file = "intePareto.CERNO.A79V.ECC.csv") 
