#GSE138695 RNA-Seq DGE Analysis pipeline 
#Valdes et. al 2023 Molecular Psychiatry Submission
#Load packages
library("edgeR")
library("limma")
library("Glimma")
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
library("data.table")
library("enrichR")
library("msigdbr")
library("RColorBrewer")
library("DESeq2")
library("topconfects")
library("BiocParallel")
library("tmod")
# Set base directory
base_dir <- getwd()
# Get ENSEMBL mart
martGRCh38.99 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                   dataset = "hsapiens_gene_ensembl",
                                   host = 'jan2020.archive.ensembl.org',
                                   path="/biomart/martservice")
GRCh38.t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id_version",
                                            "ensembl_gene_id",
                                            "hgnc_symbol"),
                             mart = martGRCh38.99)
GRCh38.t2g <- dplyr::rename(GRCh38.t2g,
                            TXNAME = ensembl_transcript_id_version,
                            ENSEMBL = ensembl_gene_id,
                            SYMBOL = hgnc_symbol)
head(GRCh38.t2g)
# set GSE138695 kallisto dir
kallisto_dir <- "~/GSE138695/kallistoOut"
# Get the samples from a txt file containing the run_accession, sample, and condition
samples <- read.table(file.path(base_dir, "GSE13895_Samples.txt"), header = TRUE, stringsAsFactors=FALSE)
# For Kallisto, describe the path to find the quant.sf files
files <- file.path(kallisto_dir, samples$SampleID, "abundance.h5")
# Apply the sample names to "files"
names(files) <- paste0(c(samples$SampleID))
# Check if all files exist
all(file.exists(files))
# Import the abundance/counts measurements using tximport
txi_lsTPM <- tximport(files, 
                      type = "kallisto", 
                      tx2gene = GRCh38.t2g, 
                      countsFromAbundance = "lengthScaledTPM")
#
sampleTable <- data.frame(SampleID=factor(rep(c(samples$SampleID))),
                          Condition=factor(rep(c(samples$condition))),
                          Subject=factor(rep(c(samples$patient))),
                          SampleName=factor(rep(c(samples$sample))))
## PERFORM DIFFERENTIAL EXPRESSION WITH LIMMA-VOOM ##
# Create DGEList
y <- DGEList(txi_lsTPM$counts,
             lib.size = colSums(txi_lsTPM$counts),
             norm.factors = calcNormFactors(txi_lsTPM$counts),
             samples = samples$sample,
             group = samples$condition)
# Add the ENTREZID, SYMBOL for each ENSEMBL ID
library(org.Hs.eg.db)
#Create a Homo Sapiens annotation
Hs_ann <- AnnotationDbi::select(org.Hs.eg.db,
                                keys=rownames(y$counts),
                                column=c("ENTREZID","SYMBOL"),
                                keytype="ENSEMBL",
                                multiVals="first")
# Remove duplicated terms
Hs_ann <- Hs_ann[!duplicated(Hs_ann[,1]),]
head(Hs_ann)
# Apply the annotation to your limma object "y"
y$genes <- Hs_ann
# View the library size for each sample
y$samples
# Load a nice color palette of 50 colors to be used for plots
myPalette <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
# convert counts to cpm and log
unfilteredExpr <- cpm(y, log=T)
# Plot the density of unfiltered gene expression for all samples within groups
par(mfrow=c(1,1))
plotDensities(unfilteredExpr, group=samples$Condition, col=myPalette[1:4], legend = "topright")
boxplot(unfilteredExpr, las=2, group=sampleTable$Disease, col=myPalette[c(1,1,1,2,2,2,3,3,3,4,4,4)], main="")
# Filter lowly expressed genes via edgeR
keep = filterByExpr(y)
table(keep)
y <- y[keep,]
# Calculate normalization factors
y <- calcNormFactors(y, method = "TMM")
# Create the group and design - rename your conditions!
Condition <- factor(sampleTable$Condition,levels=c("AD","Iso"))
Subject <- factor(sampleTable$Subject,levels=c("1","2"))
Contrasts <- paste(sampleTable$Condition, sampleTable$Subject, sep=".")
# Factorial contrast design
Contrasts <- factor(Contrasts)
Contrasts
sampleTable$Contrasts <- Contrasts
# Create the design model matrix. If you have different batches, you want want to define batches and use ~0+group+batch
# Define the model matrix 
design.C <- model.matrix(~0+Contrasts)
colnames(design.C)[c(1:4)] <- c(levels(Contrasts))
design <- model.matrix(~0+Condition+Subject)
# Create the contrast matrix for DE
contr.matrix.C <- makeContrasts(ADvIso.C = (AD.1 - Iso.1)+(AD.2 - Iso.2),
                              levels = design.C)
# See the contrast matrix
contr.matrix.C
# Can also run voomwithQualityWeights (see limma vignette) - unhash the next line if so
v.C <- voomWithQualityWeights(y, design=design.C, plot=TRUE)
# fit the linear model
fit.C <- lmFit(v.C, design.C)
# apply your contrasts
cfit.C <- contrasts.fit(fit.C, contrasts=contr.matrix.C)
# eBayes method
efit.C <- eBayes(cfit.C)
#
plotSA(efit.C, main="Final model: Mean-variance trend")
# See how many genes are differentially expressed
summary(decideTests(efit.C))
# topTable of results
res_IsoAD <- topTable(efit.C,coef = 1, adjust.method = "BH", sort.by = "t", n = Inf)
head(res_IsoAD)
write.table(res_IsoAD, file = "GSE138695.IsoAD.tsv", sep = "/t")

# CERNO enrichment analysis
library("tmod")
#Import the latest MSigDB XML File
msig <- tmodImportMSigDB("/Users/andrewbcaldwell/Library/Mobile Documents/com~apple~CloudDocs/gmtfiles/msigdb_v2023.1.Hs.xml")
#Select the Hallmark gene sets (you can selectwhichever geneset db you would like)
unique(as.character(msig$gs$Subcategory))
unique(as.character(msig$gs$Category))
Hallmark.sel <- msig$gs$Category == "H"
GOBP.sel <- msig$gs$Subcategory == "GO:BP"
GTRD.sel <- msig$gs$Subcategory == "TFT:GTRD"
KEGG.sel <- msig$gs$Subcategory == "CP:KEGG"
Wiki.sel <- msig$gs$Subcategory == "CP:WIKIPATHWAYS"
BioCARTA.sel <- msig$gs$Subcategory == "C2_CP:BIOCARTA"
Reactome.sel <- msig$gs$Subcategory == "CP:REACTOME"
miR.sel <- msig$gs$Subcategory == "MIR:MIRDB"
GOBP.exp.sel <- msig$gs$Subcategory == "C5_GO:BP"
# Load ENCODE-ChEA consensus
ECC <- read_table("~/gmtfiles/ECC_TFchipenrich.txt")
ECC <- ECC[,c(2,1)]
colnames(ECC) <- c("gene_id","go_id")
m2g_ECC <- split(ECC$gene_id, ECC$go_id)
#gt <- toTable(GOTERM)
m_ECC <- data.frame(ID=names(m2g_ECC))
m_ECC$Title <- m_ECC$ID
goset_ECC <- makeTmod(modules=m_ECC, modules2genes=m2g_ECC)
# ECC
CERNO.ECC <- tmodLimmaTest(efit.C,v.C$genes$ENTREZID, sort.by = "msd", tmodFunc = tmodCERNOtest, mset=goset_ECC)
head(CERNO.ECC)
# Hallmark
CERNO.H <- tmodLimmaTest(efit.C,v.C$genes$SYMBOL, sort.by = "msd", tmodFunc = tmodCERNOtest, mset=msig[Hallmark.sel])
CERNO.H
# Reactome
CERNO.R <- tmodLimmaTest(efit.C,v.C$genes$SYMBOL, sort.by = "msd", tmodFunc = tmodCERNOtest, mset=msig[Reactome.sel])
head(CERNO.R)
#
ECC.gmt <- read.gmt(file.path("~/gmtfiles", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt"))
ReMap.gmt <- read.gmt(file.path("~/gmtfiles", "ReMap_ChIP-seq.gmt"))
# create limma ranks (t-value) for fgsea
reslimma <- topTable(efit.C, coef=1, adjust.method="BH", sort.by = "t", n = Inf)
#Convert your limma topTable results to a data table format
reslimma_dt <- data.table(reslimma)
# Remove rows with an empty or duplicated Symbol
reslimma_dt <- subset(reslimma_dt, SYMBOL != "" )
reslimma_dt <- subset(reslimma_dt, ! duplicated(SYMBOL))
# Tidy up and order by t metric
ranks_limma <- as_tibble(reslimma_dt[order(t), list(SYMBOL, t)])
ranks_limma
# Deframe the data table
ranks <- deframe(ranks_limma)
#see the top 20 genes by t
head(ranks, 20)
#You alread imported the Gene Ontology - Biological Process as Hs.GOBP
Hs.GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
# We need to select Symbol now instead of ENTREZ ID
Hs.GOBP.Symbol <- Hs.GOBP %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.GOBP.Symbol %>% 
  head() %>% 
  lapply(head)
#Run fGSEA with with multilevel approach
#Set nproc to the (total processors available)-2
fgseaRes.GOBP <- fgsea::fgseaMultilevel(pathways=Hs.GOBP.Symbol, stats=ranks, eps = 0, nproc = 6)
# Make the Results tidy
fgseaResTidy <- fgseaRes.GOBP %>%
  as_tibble() %>%
  arrange(desc(-NES))
fgseaResTidy
#
fgseaRes.Reactome <- fgsea::fgseaMultilevel(pathways=Hs.Reactome.Symbol, stats=ranks, eps = 0, nproc = 6)
# Make the Results tidy
fgseaResTidy.R <- fgseaRes.Reactome %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy.R
#
Hs.Hallmark <- msigdbr(species = "Homo sapiens", category = "H")
Hs.Hallmark.Symbol <- Hs.Hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
fgseaRes.Hallmark <- fgsea::fgseaMultilevel(pathways=Hs.Hallmark.Symbol, stats=ranks, eps = 0, nproc = 6)
# Make the Results tidy
fgseaResTidy.H <- fgseaRes.Hallmark %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy.H
#
fgseaRes.ECC <- fgsea::fgseaMultilevel(pathways=ECC.gmt, stats=ranks, eps = 0, nproc = 6)
# Make the Results tidy
fgseaResTidy.ECC <- fgseaRes.ECC %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy.ECC

# Run Dorothea
data(dorothea_hs, package = "dorothea")
viper_regulon <- dorothea_hs %>%
  filter(confidence %in% c("A", "B","C"))
viper_regulon <- df2regulon(viper_regulon)
# Explore the regulon object
names(viper_regulon)[1:10]
viper_regulon[[1]]
#
DEsignature <- reslimma_dt
# Estimate z-score values for the GES. Check VIPER manual for details
myStatistics <- matrix(DEsignature$logFC, dimnames = list(DEsignature$SYMBOL, 'logFC') )
myPvalue <- matrix(DEsignature$adj.P.Val, dimnames = list(DEsignature$SYMBOL, 'adj.P.Val') )
mySignature <- (qnorm(myPvalue/2, lower.tail = FALSE) * sign(myStatistics))[, 1]
# Reorder and rename
generanking <- mySignature[order(mySignature, decreasing = T)]
# Estimate TF activities
mrs <- viper::msviper(ges = generanking, regulon = viper_regulon, minsize = 4, ges.filter = F, verbose = TRUE)
# See the results
summary(mrs)
# Save the results as a dataframe
TF_activities <- data.frame(Regulon = names(mrs$es$nes),
                            Size = mrs$es$size[ names(mrs$es$nes) ], 
                            NES = mrs$es$nes, 
                            p.value = mrs$es$p.value, 
                            FDR = p.adjust(mrs$es$p.value, method = 'fdr'))
# Tidy up and order by FDR
TF_activitiesTidy <- TF_activities %>%
  as_tibble() %>%
  arrange(FDR)
# Show in a nice table
DT::datatable(TF_activitiesTidy)
