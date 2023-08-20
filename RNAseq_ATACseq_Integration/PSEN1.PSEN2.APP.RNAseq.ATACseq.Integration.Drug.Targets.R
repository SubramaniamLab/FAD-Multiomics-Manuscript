#Intersection of dATAC peaks with RNA-seq data using CERNO with tmod R package
#Drug Target Analysis
#Valdes et. al 2023 Molecular Psychiatry Submission

##############################
#1. Install and Load Packages
##############################

#Install and load the following packages
BiocManager::install("DiffBind")
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("rtracklayer")
BiocManager::install("plyranges")
BiocManager::install("tmod")
BiocManager::install("GOfuncR")
BiocManager::install("qusage")
BiocManager::install("topconfects")
BiocManager::install("enrichR")
BiocManager::install("networkD3")
BiocManager::install("circlize")
BiocManager::install("svglite")

#Load the packages
library("DiffBind")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("rtracklayer")
library("plyranges")
library(GenomicFeatures)
library(dplyr)
library(tmod) 
library(GOfuncR) 
library("msigdbr")
library("qusage") 
library(topconfects)
library(enrichR)
library(nVennR)
library(networkD3)
library(circlize)
library(svglite)

###########################
#2. Import msigdb Genesets
###########################

##===========================================================================================================
#Import the latest MSigDB XML File -> download 7.1 version
msig <- tmodImportMSigDB(file.path("/Users/phoebevaldes/Desktop/ATAC-seq_2021/msigdb", "msigdb_v7.1.xml"))

#Select the GO BP gene sets (you can select whichever geneset db you would like)
Hallmark.sel <- msig$MODULES$Category == "H"

#All values come out as False
GOBP.sel <- msig$MODULES$Subcategory == "BP"

#Download msigdb genesets
msigdbr_show_species()

# OPTION #1: Use msigdbr and select which geneset you want - here is Gene Ontology - BP
Hs.GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
Hs.Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
Hs.Hallmark <- msigdbr(species = "Homo sapiens", category = "H")

#Create new ENCODE/ChEA object here
ECC.gmt <-  qusage::read.gmt(file.path("/Users/phoebevaldes/Desktop/Winter2021/LimmaVoom/gmtfiles", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt"))

#Create microarray object here
Hs.miR <- qusage::read.gmt(file.path("/Users/phoebevaldes/Desktop/Winter2021/LimmaVoom/gmtfiles", "miRTarBase_2017.gmt"))

#Make ENCODE/ChEA Consensus mset for use with tmod
ECC.MODULES <- as.data.frame(read.csv(file.path("/Users/phoebevaldes/Desktop/Winter2021/LimmaVoom/gmtfiles", "ECCTFlist_MODULES_FINAL2.csv")))
ECC.MODULES2GENES <- as.list(ECC.gmt)
rownames(ECC.MODULES) <- ECC.MODULES$ID
ECC.GENES <- as.data.frame(v2$genes$SYMBOL)
names(ECC.GENES)[1] <- "ID"
data(tmod)
msetECC <- tmod
msetECC$GENES2MODULES <- NULL
msetECC$GENES <- ECC.GENES
msetECC$MODULES2GENES <- ECC.MODULES2GENES
msetECC$MODULES <- ECC.MODULES

#Make ReMap ChIP-Seq mset for use with tmod
ReMap.gmt <- qusage::read.gmt(file.path("/Users/phoebevaldes/Desktop/Winter2021/LimmaVoom/gmtfiles", "ReMap_ChIP-seq.gmt"))
ReMap.MODULES <- as.data.frame(read.csv(file.path("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Spring2020/Rotation/gmtfiles", "ReMapTFlist_MODULES.csv")))
ReMap.MODULES2GENES <- as.list(ReMap.gmt)
rownames(ReMap.MODULES) <- ReMap.MODULES$ID
ReMap.GENES <- as.data.frame(v2$genes$SYMBOL)
names(ReMap.GENES)[1] <- "ID"
msetReMap <- tmod
msetReMap$GENES2MODULES <- NULL
msetReMap$GENES <- ReMap.GENES
msetReMap$MODULES2GENES <- ReMap.MODULES2GENES
msetReMap$MODULES <- ReMap.MODULES

#Perform wget from https://urldefense.com/v3/__https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqeoeY17M$  in TSCC
#Then perform gunzip hg38.ensGene.gtf.gz in TSCC 
base_dir <- "/Users/phoebevaldes/Desktop/ATAC-seq_2021"
txdb.hg38ens <- makeTxDbFromGFF(file.path(base_dir, "hg38.ensGene.gtf"), format = "gtf")
saveDb(txdb.hg38ens,file = "txdb.hg38ens.sqlite")
txdb.hg38ens <- loadDb("txdb.hg38ens.sqlite")
txdb.hg38ens
keytypes(txdb.hg38ens)
seqlevels(txdb.hg38ens)

##############################
#3. Get DEGs and efit object
##############################

#Calculate normalization factors
y <- calcNormFactors(y, method = "TMM")
gsva <- y

#Run voom
v <- voom(gsva, design=design, plot=TRUE)

# fit the linear model
fit <- lmFit(v)
names(fit)

# apply your contrasts
cfit <- contrasts.fit(fit, contrasts=contr.matrix)

# eBayes method
efit <- eBayes(cfit)

#Save the 'efit' object (From the previous environment saved from CEMiTool_GSVA_Networks_v5.RData file)
saveRDS(file="/Users/phoebevaldes/Desktop/ATAC-seq_2021/RDS_Files/efit.rds", efit)
efit <- readRDS(file="/Users/phoebevaldes/Desktop/ATAC-seq_2021/RDS_Files/efit.rds")

dim(efit)

#See DEGs by eBayes method
summary(decideTests(efit))


#####################################################    
#4. APP-V717I Mutation - EBAYES METHOD WITH topTABLE
#Integration with RNA-seq and ATAC-seq
#####################################################

##===========================================================================================================
##APP-V717I Mutation - EBAYES METHOD WITH topTABLE 

load("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/RData/102621_ATAC_tmod_Env.RData")

#Using the eBayes method with topTable for differential gene expression
V717I.RNA_3 <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf) #for 10.27.21 results

dim(V717I.RNA_3)
head(V717I.RNA_3)

#Use table to make volcano plots
#Write the data matrix into a text file
write.csv(V717I.RNA_3, '/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/topTable_Results/230731.topTable_results_V717Ivs.NDC.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

# Isolate the key attributes and rename the columns
sigGenes.V717I <- V717I.RNA_3 %>%
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
#2151 differential peaks
dATACpeaks.V717I.DEG 

# Isolate the key statistics of each
dATACpeaks.V717I.DEG <- dATACpeaks.V717I.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.V717I.DEG
write.csv(dATACpeaks.V717I.DEG, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATACpeaks.V717I.DEG.csv")

####################################################
#5. Do drug target analysis for APP-V717I mutation 
#for integrated dATAC peaks and DEGs - all regions
####################################################

#Read the genes file associated with drugs from Fang et al., 2022
drugTargets_Fang2022 <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Literature_Sources/Fang.etal.2022_SupplementaryTable4_DTIs.csv")

#Merge integrated dATAC peaks and DEGs with drug targets by gene symbol
drugTargets_dATACpeaks.V717I.DEG <- merge(dATACpeaks.V717I.DEG, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.V717I.DEG, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/V717I/230802.dATACpeaks.V717I.DEG.drugTargets.csv")


############################################################    
#6. APP-V717I Mutation - EBAYES METHOD WITH topTABLE
#Integration with RNA-seq and ATAC-seq - Promoter Regions
############################################################

# Split into promoter and non-promoter regions and find only unique peak intervals
dATACpeaks.V717I.DEG.promoter <- dATACpeaks.V717I.DEG[dATACpeaks.V717I.DEG$annotation == "Promoter"]
#1360 differential peaks in the promoter region
dATACpeaks.V717I.DEG.promoter
dATACpeaks.V717I.DEG.nonpromoter <- dATACpeaks.V717I.DEG[dATACpeaks.V717I.DEG$annotation != "Promoter"]
#791 differential peaks in the nonpromoter region
dATACpeaks.V717I.DEG.nonpromoter

#Write the .csv files
write.csv(dATACpeaks.V717I.DEG.promoter, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATACpeaks.V717I.DEG.promoter.csv", row.names = TRUE)
write.csv(dATACpeaks.V717I.DEG.nonpromoter, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATACpeaks.V717I.DEG.nonpromoter.csv")


# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.V717I.DEG.promoter <- unique(dATACpeaks.V717I.DEG.promoter)

V717I.dATAC.promoter.DEG.up <- dATACpeaks.V717I.DEG.promoter %>%
  dplyr::filter(da_log2FC > 0) %>%
  dplyr::filter(de_log2FC > 0)
V717I.dATAC.promoter.DEG.down <- dATACpeaks.V717I.DEG.promoter %>%
  dplyr::filter(da_log2FC < 0) %>%
  dplyr::filter(de_log2FC < 0)

#Select the increased and decreased promoter accessibility for the APP-V717I mutation
#244 peaks are differential in the enhancer (increased accesibility within the promoter region)
V717I.dATAC.promoter.DEG.up
#858 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
V717I.dATAC.promoter.DEG.down

#Write the .csv files
write.csv(V717I.dATAC.promoter.DEG.up, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATAC.promoter.DEG.up.V717I.csv")
write.csv(V717I.dATAC.promoter.DEG.down, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATAC.promoter.DEG.down.V717I.csv")

#####################################################################
#7. Do drug target analysis for APP-V717I mutation 
#for integrated dATAC peaks and DEGs - promoter up and down regions
#####################################################################

#Merge integrated dATAC peaks and DEGs with drug targets in promoter up regions
#by gene symbol
drugTargets_dATACpeaks.promoter.V717I.DEG.up <- merge(V717I.dATAC.promoter.DEG.up, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.promoter.V717I.DEG.up, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/V717I/230802.dATACpeaks.V717I.DEG.promoter.up.drugTargets.csv")

#Read .csv file with unique gene symbols with a drug target
#drugTargets_V717I.promoter.up_UNIQUE <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/V717I/230802.dATACpeaks.promoter.V717I.DEG.up.drugTargets_UNIQUE2.csv")

#Merge integrated dATAC peaks and DEGs with drug targets in promoter down regions
#by gene symbol
drugTargets_dATACpeaks.promoter.V717I.DEG.down <- merge(V717I.dATAC.promoter.DEG.down, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.promoter.V717I.DEG.down, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/V717I/230803.dATACpeaks.V717I.DEG.promoter.down.drugTargets.csv")

####################################################################    
#8. APP-V717I Mutation - EBAYES METHOD WITH topTABLE
#Integration with RNA-seq and ATAC-seq - PEREGRINE Enhancer Regions
####################################################################

# Find PEREGINE enhancer associated peaks
enhancer.PEREGRINE <- as.data.frame(read.table("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Enhancer_Data/PEREGRINE.brain.enhancer.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(enhancer.PEREGRINE)[1] <- "seqnames"
names(enhancer.PEREGRINE)[2] <- "start"
names(enhancer.PEREGRINE)[3] <- "end"

# Convert to 3 column bed
enhancer.PEREGRINE.ranges <- plyranges::as_granges(enhancer.PEREGRINE[,c(1:3)])

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
dATACpeaks.V717I.DEG.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.V717I.DEG.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
dATACpeaks.V717I.DEG.enhancer.PEREGRINE <- unique(dATACpeaks.V717I.DEG.enhancer.PEREGRINE)
V717I.dATAC.enhancer.PEREGRINE.DEG.ranges <- plyranges::reduce_ranges(dATACpeaks.V717I.DEG.enhancer.PEREGRINE)
V717I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- as.data.frame(V717I.dATAC.enhancer.PEREGRINE.DEG.ranges)
V717I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- V717I.dATAC.enhancer.PEREGRINE.DEG.ranges.df[,c(1:3)]

#Increased accessibility
V717I.dATAC.enhancer.PEREGRINE.DEG.up <- dATACpeaks.V717I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC > 0)
#Decreased accessibility
V717I.dATAC.enhancer.PEREGRINE.DEG.down <- dATACpeaks.V717I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC < 0)

#Write the .csv file
write.csv(dATACpeaks.V717I.DEG.enhancer.PEREGRINE, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATAC.enhancer.PEREGRINE.DEG.V717I.csv")

#Select the increased and decreased enhancer accessibility for the APP-V717I mutation
#133 peaks are differential in the enhancer (increased accesibility within the enhancer region)
V717I.dATAC.enhancer.PEREGRINE.DEG.up
#147 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
V717I.dATAC.enhancer.PEREGRINE.DEG.down

#Write the .csv file
write.csv(V717I.dATAC.enhancer.PEREGRINE.DEG.up, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATAC.enhancer.PEREGRINE.DEG.up.V717I.csv")
write.csv(V717I.dATAC.enhancer.PEREGRINE.DEG.down, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/V717I/230731.dATAC.enhancer.PEREGRINE.DEG.down.V717I.csv")


##################################################################
#9. Do drug target analysis for APP-V717I mutation 
#for integrated dATAC peaks and DEGs - PEREGRINE enhancer regions
##################################################################

#Merge integrated dATAC peaks and DEGs with drug targets in PEREGRINE enhancer up regions
#by gene symbol
drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.up <- merge(V717I.dATAC.enhancer.PEREGRINE.DEG.up, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.up, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/V717I/230803.dATACpeaks.V717I.DEG.enhancer.PEREGRINE.up.drugTargets.csv")


#Merge integrated dATAC peaks and DEGs with drug targets in PEREGRINE enhancer down regions
#by gene symbol
drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.down <- merge(V717I.dATAC.enhancer.PEREGRINE.DEG.down, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.down, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/V717I/230803.dATACpeaks.V717I.DEG.enhancer.PEREGRINE.down.drugTargets.csv")


#########################################################
#10. Do ranking of intePareto Results for APP-V717I mutation
#########################################################


#Read in the .csv file with the ranked gene list from intePareto for APP mutation
rankedGeneList_intePareto_V717I <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/intePareto_Results/211213.V717I.intePareto.pATAC.peregrineATAC.results.csv")

#Filter results
rankedGeneList_intePareto_V717I.df<- as.data.frame(rankedGeneList_intePareto_V717I)
rankedGeneList_intePareto_V717I.dt = subset(rankedGeneList_intePareto_V717I.df, SYMBOL != "")
rankedGeneList_intePareto_V717I.dt = subset(rankedGeneList_intePareto_V717I.df, ! duplicated(SYMBOL))

#Tidy up and order by t metric
limma_V717I <- V717I.RNA_3 #Rename variable
ranks_intePareto_V717I <- rankedGeneList_intePareto_V717I.dt

#Calculates the (nonsymmetric) set difference between gene symbols of the 
#topTable and intePareto ranks
row.names.remove.V717I <- setdiff(ranks_intePareto_V717I$SYMBOL, limma_V717I$SYMBOL) 
ranks_intePareto_V717I <- ranks_intePareto_V717I[!(ranks_intePareto_V717I$SYMBOL
                                                   %in% row.names.remove.V717I), ]
ranks_intePareto_V717I <- subset(ranks_intePareto_V717I, ! duplicated(ranks_intePareto_V717I$SYMBOL))

#####################################################################
#11. Do RNA-seq and ATAC-seq integration for APP-V717I vs. NDC 
#tmod CERNO test
#####################################################################

#GOBP - ordered intePareto Results
CERNO.GOBP.intePareto_V717I <- tmodCERNOtest(ranks_intePareto_V717I$SYMBOL, 
                                             modules = NULL,
                                             qval = 0.05,
                                             order.by = "pval",
                                             filter = FALSE,
                                             mset = msig[GOBP.sel],
                                             cols = "Title")

#Write CERNO GOBP results to .csv file
write.csv(CERNO.GOBP.intePareto_V717I, file = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/V717I/230802.V717I.dATAC.intePareto.DEG.tmodCERNO.GOBP.csv")

#GOBP - ordered intePareto Results
CERNO.Hallmark.intePareto_V717I <- tmodCERNOtest(ranks_intePareto_V717I$SYMBOL, 
                                             modules = NULL,
                                             qval = 0.05,
                                             order.by = "pval",
                                             filter = FALSE,
                                             mset = msig[Hallmark.sel],
                                             cols = "Title")

#Write CERNO Hallmark results to .csv file
write.csv(CERNO.Hallmark.intePareto_V717I, file = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/V717I/230802.V717I.dATAC.intePareto.DEG.tmodCERNO.Hallmark.csv")


########################################################
#12. Do drug target analysis for APP-V717I mutation 
#for CERNO results from intePareto Results
########################################################

#Read drug targets by CARDO class from Cummings et al., 2021
CADROclasses.Cummings2021 <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Literature_Sources/Cummings_etal.,2022_Drug_Agents_CADRO_Classes.csv")

#Read APP-V717I CERNO results curated by Common Alzheimer's Disease Research Ontology (CADRO)
#mechanism classes from Cummings et al., 2021
CERNO_CADROclasses.V717I <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/V717I/211214_intepareto_CERNO_all_V717I.csv")

#Merge by CADRO mechanism class
CERNO_CADROclasses.drugTargets.V717I <- merge(CERNO_CADROclasses.V717I, CADROclasses.Cummings2021, by = "CADRO_mechanism_class")

#Write the .csv file
write.csv(CERNO_CADROclasses.drugTargets.V717I, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_Pathway_Drugs_Results/V717I/230802.CERNOpathways_drugTargets.V717I.csv")

##============================================================================================================

###################################################    
#13. PSEN1-A79V Mutation - EBAYES METHOD WITH topTABLE  
###################################################

#Using the eBayes method with topTable for differential gene expression
A79V.RNA_3 <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf) #for 10.27.21 results

dim(A79V.RNA_3)
head(A79V.RNA_3)

#Use table to make volcano plots
#Write the data matrix into a text file
write.csv(A79V.RNA_3, '/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/topTable_Results/230731.topTable_results_A79Vvs.NDC.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

# Isolate the key attributes and rename the columns
sigGenes.A79V <- A79V.RNA_3 %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.A79V)

# Remove duplicated ENSEMBL terms
sigGenes.A79V = subset(sigGenes.A79V, ENSEMBL != "" )
sigGenes.A79V = subset(sigGenes.A79V, ! duplicated(ENSEMBL))
dim(sigGenes.A79V)
head(sigGenes.A79V)

# Match the genes from gene expression data to the differentially accessible ATAC peaks
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
#4253 differential peaks
dATACpeaks.A79V.DEG 

# Isolate the key statistics of each
dATACpeaks.A79V.DEG <- dATACpeaks.A79V.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.A79V.DEG
write.csv(dATACpeaks.A79V.DEG, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATACpeaks.A79V.DEG.csv")

#################################################
#14. Do drug target analysis for PSEN1-A79V mutation
#for integrated dATAC peaks and DEGs
#################################################

#Merge integrated dATAC peaks and DEGs with drug targets by gene symbol
drugTargets_dATACpeaks.A79V.DEG <- merge(dATACpeaks.A79V.DEG, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.A79V.DEG, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/A79V/230802.dATACpeaks.A79V.DEG.drugTargets.csv")


############################################################    
#15. PSEN1-A79V Mutation - EBAYES METHOD WITH topTABLE
#Integration with RNA-seq and ATAC-seq - Promoter Regions
############################################################

# Split into promoter and non-promoter regions and find only unique peak intervals
dATACpeaks.A79V.DEG.promoter <- dATACpeaks.A79V.DEG[dATACpeaks.A79V.DEG$annotation == "Promoter"]
#2740 differential peaks in the promoter region
dATACpeaks.A79V.DEG.promoter
dATACpeaks.A79V.DEG.nonpromoter <- dATACpeaks.A79V.DEG[dATACpeaks.A79V.DEG$annotation != "Promoter"]
#1513 differential peaks in the nonpromoter region
dATACpeaks.A79V.DEG.nonpromoter

#Write the .csv files
write.csv(dATACpeaks.A79V.DEG.promoter, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATACpeaks.A79V.DEG.promoter.csv", row.names = TRUE)
write.csv(dATACpeaks.A79V.DEG.nonpromoter, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATACpeaks.A79V.DEG.nonpromoter.csv")


# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.A79V.DEG.promoter <- unique(dATACpeaks.A79V.DEG.promoter)

A79V.dATAC.promoter.DEG.up <- dATACpeaks.A79V.DEG.promoter %>%
  dplyr::filter(da_log2FC > 0) %>%
  dplyr::filter(de_log2FC > 0)
A79V.dATAC.promoter.DEG.down <- dATACpeaks.A79V.DEG.promoter %>%
  dplyr::filter(da_log2FC < 0) %>%
  dplyr::filter(de_log2FC < 0)

#Select the increased and decreased promoter accessibility for the PSEN1-A79V mutation
#608 peaks are differential in the enhancer (increased accessibility within the promoter region)
A79V.dATAC.promoter.DEG.up
#1064 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
A79V.dATAC.promoter.DEG.down

#Write the .csv files
write.csv(A79V.dATAC.promoter.DEG.up, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATAC.promoter.DEG.up.A79V.csv")
write.csv(A79V.dATAC.promoter.DEG.down, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATAC.promoter.DEG.down.A79V.csv")

#####################################################################
#16. Do drug target analysis for PSEN1-A79V mutation 
#for integrated dATAC peaks and DEGs - promoter up and down regions
#####################################################################

#Merge integrated dATAC peaks and DEGs with drug targets in promoter up regions
#by gene symbol
drugTargets_dATACpeaks.promoter.A79V.DEG.up <- merge(A79V.dATAC.promoter.DEG.up, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.promoter.A79V.DEG.up, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/A79V/230803.dATACpeaks.A79V.DEG.promoter.up.drugTargets.csv")


#Merge integrated dATAC peaks and DEGs with drug targets in promoter down regions
#by gene symbol
drugTargets_dATACpeaks.promoter.A79V.DEG.down <- merge(A79V.dATAC.promoter.DEG.down, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.promoter.A79V.DEG.down, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/A79V/230803.dATACpeaks.A79V.DEG.promoter.down.drugTargets.csv")



####################################################################    
#17. PSEN1-A79V Mutation - EBAYES METHOD WITH topTABLE
#Integration with RNA-seq and ATAC-seq - PEREGRINE Enhancer Regions
####################################################################

# Find PEREGINE enhancer associated peaks
enhancer.PEREGRINE <- as.data.frame(read.table("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Enhancer_Data/PEREGRINE.brain.enhancer.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(enhancer.PEREGRINE)[1] <- "seqnames"
names(enhancer.PEREGRINE)[2] <- "start"
names(enhancer.PEREGRINE)[3] <- "end"

# Convert to 3 column bed
enhancer.PEREGRINE.ranges <- plyranges::as_granges(enhancer.PEREGRINE[,c(1:3)])

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
dATACpeaks.A79V.DEG.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.A79V.DEG.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
dATACpeaks.A79V.DEG.enhancer.PEREGRINE <- unique(dATACpeaks.A79V.DEG.enhancer.PEREGRINE)
A79V.dATAC.enhancer.PEREGRINE.DEG.ranges <- plyranges::reduce_ranges(dATACpeaks.A79V.DEG.enhancer.PEREGRINE)
A79V.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- as.data.frame(A79V.dATAC.enhancer.PEREGRINE.DEG.ranges)
A79V.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- A79V.dATAC.enhancer.PEREGRINE.DEG.ranges.df[,c(1:3)]

#Increased accessibility
A79V.dATAC.enhancer.PEREGRINE.DEG.up <- dATACpeaks.A79V.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC > 0)
#Decreased accessibility
A79V.dATAC.enhancer.PEREGRINE.DEG.down <- dATACpeaks.A79V.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC < 0)

#Write the .csv file
write.csv(dATACpeaks.A79V.DEG.enhancer.PEREGRINE, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATAC.enhancer.PEREGRINE.DEG.A79V.csv")

#Select the increased and decreased enhancer accessibility for the APP-V717I mutation
#355 peaks are differential in the enhancer (increased accesibility within the enhancer region)
A79V.dATAC.enhancer.PEREGRINE.DEG.up
#261 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
A79V.dATAC.enhancer.PEREGRINE.DEG.down

#Write the .csv file
write.csv(A79V.dATAC.enhancer.PEREGRINE.DEG.up, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATAC.enhancer.PEREGRINE.DEG.up.A79V.csv")
write.csv(A79V.dATAC.enhancer.PEREGRINE.DEG.down, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/A79V/230731.dATAC.enhancer.PEREGRINE.DEG.down.A79V.csv")


###################################################################
#18. Do drug target analysis for PSEN1-A79V mutation 
#for integrated dATAC peaks and DEGs - PEREGRINE enhancer regions
###################################################################

#Merge integrated dATAC peaks and DEGs with drug targets in PEREGRINE enhancer up regions
#by gene symbol
drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.up <- merge(A79V.dATAC.enhancer.PEREGRINE.DEG.up, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.up, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/A79V/230803.dATACpeaks.A79V.DEG.enhancer.PEREGRINE.up.drugTargets.csv")


#Merge integrated dATAC peaks and DEGs with drug targets in PEREGRINE enhancer down regions
#by gene symbol
drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.down <- merge(A79V.dATAC.enhancer.PEREGRINE.DEG.down, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.down, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/A79V/230803.dATACpeaks.A79V.DEG.enhancer.PEREGRINE.down.drugTargets.csv")


###########################################################
##19. Do ranking of intePareto Results for PSEN1-A79V Mutation
###########################################################

#Read in the .csv file with the ranked gene list from intePareto for APP mutation
rankedGeneList_intePareto_A79V <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/intePareto_Results/211213.A79V.intePareto.pATAC.peregrineATAC.results.csv")

#Filter results
rankedGeneList_intePareto_A79V.df<- as.data.frame(rankedGeneList_intePareto_A79V)
rankedGeneList_intePareto_A79V.dt = subset(rankedGeneList_intePareto_A79V.df, SYMBOL != "")
rankedGeneList_intePareto_A79V.dt = subset(rankedGeneList_intePareto_A79V.df, ! duplicated(SYMBOL))

#Tidy up and order by t metric
limma_A79V <- A79V.RNA_3 #Rename variable
ranks_intePareto_A79V <- rankedGeneList_intePareto_A79V.dt

#Calculates the (nonsymmetric) set difference between gene symbols of the 
#topTable and intePareto ranks
row.names.remove.A79V <- setdiff(ranks_intePareto_A79V$SYMBOL, limma_A79V$SYMBOL) 
ranks_intePareto_A79V <- ranks_intePareto_A79V[!(ranks_intePareto_A79V$SYMBOL
                                                   %in% row.names.remove.A79V), ]
ranks_intePareto_A79V <- subset(ranks_intePareto_A79V, ! duplicated(ranks_intePareto_A79V$SYMBOL))

#####################################################################
#20. Do RNA-seq and ATAC-seq integration for PSEN1-A79V vs. NDC 
#tmod CERNO test
#####################################################################

#GOBP - ordered intePareto Results
CERNO.GOBP.intePareto_A79V <- tmodCERNOtest(ranks_intePareto_A79V$SYMBOL, 
                                             modules = NULL,
                                             qval = 0.05,
                                             order.by = "pval",
                                             filter = FALSE,
                                             mset = msig[GOBP.sel],
                                             cols = "Title")

#Write CERNO GOBP results to .csv file
write.csv(CERNO.GOBP.intePareto_A79V, file = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/A79V/230802.A79V.dATAC.intePareto.DEG.tmodCERNO.GOBP.csv")

#GOBP - ordered intePareto Results
CERNO.Hallmark.intePareto_A79V <- tmodCERNOtest(ranks_intePareto_A79V$SYMBOL, 
                                                 modules = NULL,
                                                 qval = 0.05,
                                                 order.by = "pval",
                                                 filter = FALSE,
                                                 mset = msig[Hallmark.sel],
                                                 cols = "Title")

#Write CERNO Hallmark results to .csv file
write.csv(CERNO.Hallmark.intePareto_A79V, file = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/A79V/230802.A79V.dATAC.intePareto.DEG.tmodCERNO.Hallmark.csv")


########################################################
#21. Do drug target analysis for PSEN1-A79V mutation 
#for CERNO results from intePareto Results
########################################################

#Read PSEN1-A79V CERNO results curated by Common Alzheimer's Disease Research Ontology (CADRO)
#mechanism classes from Cummings et al., 2021
CERNO_CADROclasses.A79V <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/A79V/211214_intepareto_CERNO_all_A79V.csv")

#Merge by CADRO mechanism class
CERNO_CADROclasses.drugTargets.A79V <- merge(CERNO_CADROclasses.A79V, CADROclasses.Cummings2021, by = "CADRO_mechanism_class")

#Write the .csv file
write.csv(CERNO_CADROclasses.drugTargets.A79V, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_Pathway_Drugs_Results/A79V/230802.CERNOpathways_drugTargets.A79V.csv")

##============================================================================================================

###################################################    
#22. PSEN2-N141I Mutation - EBAYES METHOD WITH topTABLE  
###################################################

#Using the eBayes method with topTable for differential gene expression
N141I.RNA_3 <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

dim(N141I.RNA_3)
head(N141I.RNA_3)

#Use table to make volcano plots
#Write the data matrix into a text file
write.csv(N141I.RNA_3, '/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/topTable_Results/230801.topTable_results_N141Ivs.NDC.csv', append=FALSE, sep ="\t", dec = ".", row.names= TRUE, col.names = TRUE)

# Isolate the key attributes and rename the columns
sigGenes.N141I <- N141I.RNA_3 %>%
  dplyr::select(ENSEMBL, de_log2FC = logFC, de_padj = adj.P.Val)
dim(sigGenes.N141I)

# Remove duplicated ENSEMBL terms
sigGenes.N141I = subset(sigGenes.N141I, ENSEMBL != "" )
sigGenes.N141I = subset(sigGenes.N141I, ! duplicated(ENSEMBL))
dim(sigGenes.N141I)
head(sigGenes.N141I)

# Match the genes from gene expression data to the differentially accessible ATAC peaks
dATACpeaks.N141I.w <- dATACpeaks.A141I.ranges[na.omit(match(rownames(sigGenes.N141I), names(dATACpeaks.A141I.ranges)))]
# Add the differential expression statistics of genes ID's by RNA-Seq and ATAC-Seq
mcols(dATACpeaks.N141I.w) <- sigGenes.N141I[match(names(dATACpeaks.N141I.w), rownames(sigGenes.N141I)),]

# Overlap the dATAC with dE peak regions
dATACpeaks.A141I <- dATACpeaks.A141I %>%
  join_overlap_left(dATACpeaks.N141I.w)

# Filter for both differential gene expression and differential accessibility
dATACpeaks.N141I.DEG  <- dATACpeaks.A141I  %>%
  dplyr::filter(de_padj < 0.05) %>%
  dplyr::filter(abs(de_log2FC) > 0) %>%
  dplyr::filter(da_padj < 0.05) %>%
  dplyr::filter(abs(da_log2FC) > 0)

# View the final granges object of dATAC peaks with differential gene expression
#2632 differential peaks
dATACpeaks.N141I.DEG 

# Isolate the key statistics of each
dATACpeaks.N141I.DEG <- dATACpeaks.N141I.DEG[, c("annotation", "SYMBOL", "ENTREZID", "da_log2FC", "da_padj", "de_log2FC", "de_padj")] 
dATACpeaks.N141I.DEG
write.csv(dATACpeaks.N141I.DEG, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATACpeaks.N141I.DEG.csv")

##################################################
#23. Do drug target analysis for PSEN2-N141I mutation
#for integrated dATAC peaks and DEGs
##################################################

#Merge integrated dATAC peaks and DEGs with drug targets by gene symbol
drugTargets_dATACpeaks.N141I.DEG <- merge(dATACpeaks.N141I.DEG, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.N141I.DEG, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/N141I/230802.dATACpeaks.N141I.DEG.drugTargets.csv")


############################################################    
#24. PSEN2-N141I Mutation - EBAYES METHOD WITH topTABLE
#Integration with RNA-seq and ATAC-seq - Promoter Regions
############################################################

# Split into promoter and non-promoter regions and find only unique peak intervals
dATACpeaks.N141I.DEG.promoter <- dATACpeaks.N141I.DEG[dATACpeaks.N141I.DEG$annotation == "Promoter"]
#1611 differential peaks in the promoter region
dATACpeaks.N141I.DEG.promoter
dATACpeaks.N141I.DEG.nonpromoter <- dATACpeaks.N141I.DEG[dATACpeaks.N141I.DEG$annotation != "Promoter"]
#1021 differential peaks in the nonpromoter region
dATACpeaks.N141I.DEG.nonpromoter

#Write the .csv files
write.csv(dATACpeaks.N141I.DEG.promoter, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATACpeaks.N141I.DEG.promoter.csv", row.names = TRUE)
write.csv(dATACpeaks.N141I.DEG.nonpromoter, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATACpeaks.N141I.DEG.nonpromoter.csv")


# Split into increased and decreased accessibility for promoter peaks
dATACpeaks.N141I.DEG.promoter <- unique(dATACpeaks.N141I.DEG.promoter)

N141I.dATAC.promoter.DEG.up <- dATACpeaks.N141I.DEG.promoter %>%
  dplyr::filter(da_log2FC > 0) %>%
  dplyr::filter(de_log2FC > 0)
N141I.dATAC.promoter.DEG.down <- dATACpeaks.N141I.DEG.promoter %>%
  dplyr::filter(da_log2FC < 0) %>%
  dplyr::filter(de_log2FC < 0)

#Select the increased and decreased promoter accessibility for the PSEN1-A79V mutation
#298 peaks are differential in the enhancer (increased accessibility within the promoter region)
N141I.dATAC.promoter.DEG.up
#440 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
N141I.dATAC.promoter.DEG.down

#Write the .csv files
write.csv(N141I.dATAC.promoter.DEG.up, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATAC.promoter.DEG.up.N141I.csv")
write.csv(N141I.dATAC.promoter.DEG.down, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATAC.promoter.DEG.down.N141I.csv")

#####################################################################
#25. Do drug target analysis for PSEN2-N141I mutation 
#for integrated dATAC peaks and DEGs - promoter up and down regions
#####################################################################

#Merge integrated dATAC peaks and DEGs with drug targets in promoter up regions
#by gene symbol
drugTargets_dATACpeaks.promoter.N141I.DEG.up <- merge(N141I.dATAC.promoter.DEG.up, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.promoter.N141I.DEG.up, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/N141I/230803.dATACpeaks.N141I.DEG.promoter.up.drugTargets.csv")


#Merge integrated dATAC peaks and DEGs with drug targets in promoter down regions
#by gene symbol
drugTargets_dATACpeaks.promoter.N141I.DEG.down <- merge(N141I.dATAC.promoter.DEG.down, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.promoter.N141I.DEG.down, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/N141I/230803.dATACpeaks.N141I.DEG.promoter.down.drugTargets.csv")


####################################################################    
#26. PSEN2-N141I Mutation - EBAYES METHOD WITH topTABLE
#Integration with RNA-seq and ATAC-seq - PEREGRINE Enhancer Regions
####################################################################

# Find PEREGINE enhancer associated peaks
enhancer.PEREGRINE <- as.data.frame(read.table("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Enhancer_Data/PEREGRINE.brain.enhancer.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
names(enhancer.PEREGRINE)[1] <- "seqnames"
names(enhancer.PEREGRINE)[2] <- "start"
names(enhancer.PEREGRINE)[3] <- "end"

# Convert to 3 column bed
enhancer.PEREGRINE.ranges <- plyranges::as_granges(enhancer.PEREGRINE[,c(1:3)])

# Overlap dATAC nonpromoter peaks with PEREGRINE enhancer regions
dATACpeaks.N141I.DEG.enhancer.PEREGRINE <- plyranges::join_overlap_inner(dATACpeaks.N141I.DEG.nonpromoter,enhancer.PEREGRINE.ranges)

# Split into increased and decreased accessibility
dATACpeaks.N141I.DEG.enhancer.PEREGRINE <- unique(dATACpeaks.N141I.DEG.enhancer.PEREGRINE)
N141I.dATAC.enhancer.PEREGRINE.DEG.ranges <- plyranges::reduce_ranges(dATACpeaks.N141I.DEG.enhancer.PEREGRINE)
N141I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- as.data.frame(N141I.dATAC.enhancer.PEREGRINE.DEG.ranges)
N141I.dATAC.enhancer.PEREGRINE.DEG.ranges.df <- N141I.dATAC.enhancer.PEREGRINE.DEG.ranges.df[,c(1:3)]

#Increased accessibility
N141I.dATAC.enhancer.PEREGRINE.DEG.up <- dATACpeaks.N141I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC > 0)
#Decreased accessibility
N141I.dATAC.enhancer.PEREGRINE.DEG.down <- dATACpeaks.N141I.DEG.enhancer.PEREGRINE %>%
  dplyr::filter(de_log2FC < 0)

#Write the .csv file
write.csv(dATACpeaks.N141I.DEG.enhancer.PEREGRINE, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATAC.enhancer.PEREGRINE.DEG.N141I.csv")

#Select the increased and decreased enhancer accessibility for the APP-V717I mutation
#233 peaks are differential in the enhancer (increased accesibility within the enhancer region)
N141I.dATAC.enhancer.PEREGRINE.DEG.up
#137 peaks are differential in the enhancer (decreased accessibility within the enchancer reigon)
N141I.dATAC.enhancer.PEREGRINE.DEG.down

#Write the .csv file
write.csv(N141I.dATAC.enhancer.PEREGRINE.DEG.up, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATAC.enhancer.PEREGRINE.DEG.up.N141I.csv")
write.csv(N141I.dATAC.enhancer.PEREGRINE.DEG.down, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Integration_Results/N141I/230801.dATAC.enhancer.PEREGRINE.DEG.down.N141I.csv")

###################################################################
#27. Do drug target analysis for PSEN1-A79V mutation 
#for integrated dATAC peaks and DEGs - PEREGRINE enhancer regions
###################################################################

#Merge integrated dATAC peaks and DEGs with drug targets in PEREGRINE enhancer up regions
#by gene symbol
drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.up <- merge(N141I.dATAC.enhancer.PEREGRINE.DEG.up, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.up, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/N141I/230803.dATACpeaks.N141I.DEG.enhancer.PEREGRINE.up.drugTargets.csv")


#Merge integrated dATAC peaks and DEGs with drug targets in PEREGRINE enhancer down regions
#by gene symbol
drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.down <- merge(N141I.dATAC.enhancer.PEREGRINE.DEG.down, drugTargets_Fang2022, by = "SYMBOL")

#Write the .csv file
write.csv(drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.down, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_DEG_Drugs_Results/N141I/230803.dATACpeaks.N141I.DEG.enhancer.PEREGRINE.down.drugTargets.csv")


############################################################
#28. Do ranking of intePareto Results for PSEN2-N141I Mutation
############################################################

#Read in the .csv file with the ranked gene list from intePareto for APP mutation
rankedGeneList_intePareto_N141I <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/intePareto_Results/211213.N141I.intePareto.pATAC.peregrineATAC.results.csv")

#Filter results
rankedGeneList_intePareto_N141I.df<- as.data.frame(rankedGeneList_intePareto_N141I)
rankedGeneList_intePareto_N141I.dt = subset(rankedGeneList_intePareto_N141I.df, SYMBOL != "")
rankedGeneList_intePareto_N141I.dt = subset(rankedGeneList_intePareto_N141I.df, ! duplicated(SYMBOL))

#Tidy up and order by t metric
limma_N141I <- N141I.RNA_3 #Rename variable
ranks_intePareto_N141I <- rankedGeneList_intePareto_N141I.dt

#Calculates the (nonsymmetric) set difference between gene symbols of the 
#topTable and intePareto ranks
row.names.remove.N141I <- setdiff(ranks_intePareto_N141I$SYMBOL, limma_N141I$SYMBOL) 
ranks_intePareto_N141I <- ranks_intePareto_N141I[!(ranks_intePareto_N141I$SYMBOL
                                                 %in% row.names.remove.N141I), ]
ranks_intePareto_N141I <- subset(ranks_intePareto_N141I, ! duplicated(ranks_intePareto_N141I$SYMBOL))

#####################################################################
#29. Do RNA-seq and ATAC-seq integration for PSEN1-A79V vs. NDC 
#tmod CERNO test
#####################################################################

#GOBP - ordered intePareto Results
CERNO.GOBP.intePareto_N141I <- tmodCERNOtest(ranks_intePareto_N141I$SYMBOL, 
                                            modules = NULL,
                                            qval = 0.05,
                                            order.by = "pval",
                                            filter = FALSE,
                                            mset = msig[GOBP.sel],
                                            cols = "Title")

#Write CERNO GOBP results to .csv file
write.csv(CERNO.GOBP.intePareto_N141I, file = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/N141I/230802.N141I.dATAC.intePareto.DEG.tmodCERNO.GOBP.csv")

#GOBP - ordered intePareto Results
CERNO.Hallmark.intePareto_N141I <- tmodCERNOtest(ranks_intePareto_N141I$SYMBOL, 
                                                modules = NULL,
                                                qval = 0.05,
                                                order.by = "pval",
                                                filter = FALSE,
                                                mset = msig[Hallmark.sel],
                                                cols = "Title")

#Write CERNO Hallmark results to .csv file
write.csv(CERNO.Hallmark.intePareto_N141I, file = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/N141I/230802.N141I.dATAC.intePareto.DEG.tmodCERNO.Hallmark.csv")


########################################################
#30. Do drug target analysis for PSEN2-N141I mutation 
#for CERNO results from integrated dATAC peaks and DEGs
########################################################

#Read PSEN1-A79V CERNO results curated by Common Alzheimer's Disease Research Ontology (CADRO)
#mechanism classes from Cummings et al., 2021
CERNO_CADROclasses.N141I <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/tmodCERNO_Results/N141I/211214_intepareto_CERNO_all_N141I.csv")

#Merge by CADRO mechanism class
CERNO_CADROclasses.drugTargets.N141I <- merge(CERNO_CADROclasses.N141I, CADROclasses.Cummings2021, by = "CADRO_mechanism_class")

#Write the .csv file
write.csv(CERNO_CADROclasses.drugTargets.N141I, file ="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Integrated_Pathway_Drugs_Results/N141I/230802.CERNOpathways_drugTargets.N141I.csv")

##===================================================================================================================================================

###########################################################
#31. Get the common drug targets by ID for promoter up regions
###########################################################

#Get common drug targets by drugID
APP_PSEN1_common_drugTargets.promoter.up.DEGs <- merge(drugTargets_dATACpeaks.promoter.V717I.DEG.up, drugTargets_dATACpeaks.promoter.A79V.DEG.up, by = "Drug_IDs")
PSEN1_PSEN2_common_drugTargets.promoter.up.DEGs <- merge(drugTargets_dATACpeaks.promoter.A79V.DEG.up, drugTargets_dATACpeaks.promoter.N141I.DEG.up, by = "Drug_IDs")
APP_PSEN2_common_drugTargets.promoter.up.DEGs <- merge(drugTargets_dATACpeaks.promoter.V717I.DEG.up, drugTargets_dATACpeaks.promoter.N141I.DEG.up, by = "Drug_IDs")
APP_PSEN1_PSEN2_common_drugTargets.promoter.up.DEGs <- merge(APP_PSEN1_common_drugTargets.promoter.up.DEGs, drugTargets_dATACpeaks.promoter.N141I.DEG.up, by = "Drug_IDs")


#Write .csv files
write.csv(APP_PSEN1_common_drugTargets.promoter.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Up.DEGs/230803_Common_Drug_Targets_Promoter.Up_APP_PSEN1.csv", sep="\t", quote=F, row.names=T)
write.csv(PSEN1_PSEN2_common_drugTargets.promoter.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Up.DEGs/230803_Common_Drug_Targets_Promoter.Up_PSEN1_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN2_common_drugTargets.promoter.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Up.DEGs/230803_Common_Drug_Targets_Promoter.Up_APP_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN1_PSEN2_common_drugTargets.promoter.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Up.DEGs/230803_Common_Drug_Targets_Promoter.Up_allMutations.csv", sep="\t", quote=F, row.names=T)

#################################################################
#32. Create Venn Diagram for drug targets based on integrated
#DEGs and dATAC peaks - promoter up regions
#################################################################

#Create drug target list for APP-V717I 
drugTargets_V717I.promoter.up_list <- list(drugTargets_dATACpeaks.promoter.V717I.DEG.up = as.character(unique(drugTargets_dATACpeaks.promoter.V717I.DEG.up$Drug_IDs)))

#Create drug target list for PSEN1-A79V 
drugTargets_A79V.promoter.up_list <- list(drugTargets_dATACpeaks.promoter.A79V.DEG.up = as.character(unique(drugTargets_dATACpeaks.promoter.A79V.DEG.up$Drug_IDs)))

#Create drug target list for PSEN2-N141I 
drugTargets_N141I.promoter.up_list <- list(drugTargets_dATACpeaks.promoter.N141I.DEG.up = as.character(unique(drugTargets_dATACpeaks.promoter.N141I.DEG.up$Drug_IDs)))

#Plot the nVenn Diagram 
myNV_dATAC_DEGs_drugTargets_promoter.up <- plotVenn(list(drugTargets_V717I.promoter.up_list, drugTargets_A79V.promoter.up_list, drugTargets_N141I.promoter.up_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)

#Show the .svg file
showSVG(myNV_dATAC_DEGs_drugTargets_promoter.up, opacity=0.3,outFile = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_promoter.up.DEGs_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))

######################################################################
#33. Create Chord Diagram for subset of drug targets across all mutations
#Promoter Up Regions
######################################################################

#Get unique gene symbols across common drug targets of all mutations
#to make input matrix with 0's and 1's
chosenDrugTargets.Promoter.Up <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Up.DEGs/230803_Chosen_Drug_Targets_Promoter.Up_allMutations.csv")

#Filter by drug and get unique symbols
uniqueGenes.DB00398.Promoter.Up1 <- chosenDrugTargets.Promoter.Up %>%
  dplyr::filter(Drug_IDs == "DB00398") 
uniqueGenes.DB00398.Promoter.Up.1a <- data.frame(unique(uniqueGenes.DB00398.Promoter.Up1$SYMBOL.x))
uniqueGenes.DB00398.Promoter.Up.1b <- data.frame(unique(uniqueGenes.DB00398.Promoter.Up1$SYMBOL.y))
uniqueGenes.DB00398.Promoter.Up.1c <- data.frame(unique(uniqueGenes.DB00398.Promoter.Up1$SYMBOL))

uniqueGenes.DB01254.Promoter.Up2 <- chosenDrugTargets.Promoter.Up %>%
  dplyr::filter(Drug_IDs == "DB01254") 
uniqueGenes.DB01254.Promoter.Up.2a <- data.frame(unique(uniqueGenes.DB01254.Promoter.Up2$SYMBOL.x))
uniqueGenes.DB01254.Promoter.Up.2b <- data.frame(unique(uniqueGenes.DB01254.Promoter.Up2$SYMBOL.y))
uniqueGenes.DB01254.Promoter.Up.2c <- data.frame(unique(uniqueGenes.DB01254.Promoter.Up2$SYMBOL))

uniqueGenes.DB01268.Promoter.Up3 <- chosenDrugTargets.Promoter.Up %>%
  dplyr::filter(Drug_IDs == "DB01268") 
uniqueGenes.DB01268.Promoter.Up.3a <- data.frame(unique(uniqueGenes.DB01268.Promoter.Up3$SYMBOL.x))
uniqueGenes.DB01268.Promoter.Up.3b <- data.frame(unique(uniqueGenes.DB01268.Promoter.Up3$SYMBOL.y))
uniqueGenes.DB01268.Promoter.Up.3c <- data.frame(unique(uniqueGenes.DB01268.Promoter.Up3$SYMBOL))


uniqueGenes.DB06616.Promoter.Up4 <- chosenDrugTargets.Promoter.Up %>%
  dplyr::filter(Drug_IDs == "DB06616") 
uniqueGenes.DB06616.Promoter.Up.4a <- data.frame(unique(uniqueGenes.DB06616.Promoter.Up4$SYMBOL.x))
uniqueGenes.DB06616.Promoter.Up.4b <- data.frame(unique(uniqueGenes.DB06616.Promoter.Up4$SYMBOL.y))
uniqueGenes.DB06616.Promoter.Up.4c <- data.frame(unique(uniqueGenes.DB06616.Promoter.Up4$SYMBOL))


uniqueGenes.DB06626.Promoter.Up5 <- chosenDrugTargets.Promoter.Up %>%
  dplyr::filter(Drug_IDs == "DB06626") 
uniqueGenes.DB06626.Promoter.Up.5a <- data.frame(unique(uniqueGenes.DB06626.Promoter.Up5$SYMBOL.x))
uniqueGenes.DB06626.Promoter.Up.5b <- data.frame(unique(uniqueGenes.DB06626.Promoter.Up5$SYMBOL.y))
uniqueGenes.DB06626.Promoter.Up.5c <- data.frame(unique(uniqueGenes.DB06626.Promoter.Up5$SYMBOL))

uniqueGenes.DB09079.Promoter.Up6 <- chosenDrugTargets.Promoter.Up %>%
  dplyr::filter(Drug_IDs == "DB09079") 
uniqueGenes.DB09079.Promoter.Up.6a <- data.frame(unique(uniqueGenes.DB09079.Promoter.Up6$SYMBOL.x))
uniqueGenes.DB09079.Promoter.Up.6b <- data.frame(unique(uniqueGenes.DB09079.Promoter.Up6$SYMBOL.y))
uniqueGenes.DB09079.Promoter.Up.6c <- data.frame(unique(uniqueGenes.DB09079.Promoter.Up6$SYMBOL))

#Read .csv file with drug values across PSEN1/PSEN2/APP
chordDiagram_Values_Promoter.Up <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Up.DEGs/230803_Promoter.Up_chordDiag.Values_allMutations.csv")

#Source: https://urldefense.com/v3/__https://stackoverflow.com/questions/73609132/how-to-create-a-chord-diagram-in-r__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqcD_PAlf$ 
#Set the colors for each drug
#Source: https://urldefense.com/v3/__https://r-charts.com/flow/chord-diagram/__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqe_Um70Q$ 
colors <- c(Sorafenib = "#CC79A7", Dasatinib = "#009e73",
            Sunitinib = "#449AE4", Bosutinib = "#871C9A",
            Axitinib = "#FD8305", Nintedanib = "#E64A00")

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Promoter.Up_Chord_labeled.svg")

set.seed(123)
chordDiagram(chordDiagram_Values_Promoter.Up, grid.col = colors)
dev.off()

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Promoter.Up_Chord_unlabeled.svg")

set.seed(123)
chordDiagram(chordDiagram_Values_Promoter.Up, grid.col = colors, annotationTrack = "grid")
dev.off()


###########################################################
#34. Get the common drug targets by ID for promoter down regions
###########################################################

#Get common drug targets by drugID
APP_PSEN1_common_drugTargets.promoter.down.DEGs <- merge(drugTargets_dATACpeaks.promoter.V717I.DEG.down, drugTargets_dATACpeaks.promoter.A79V.DEG.down, by = "Drug_IDs")
PSEN1_PSEN2_common_drugTargets.promoter.down.DEGs <- merge(drugTargets_dATACpeaks.promoter.A79V.DEG.down, drugTargets_dATACpeaks.promoter.N141I.DEG.down, by = "Drug_IDs")
APP_PSEN2_common_drugTargets.promoter.down.DEGs <- merge(drugTargets_dATACpeaks.promoter.V717I.DEG.down, drugTargets_dATACpeaks.promoter.N141I.DEG.down, by = "Drug_IDs")
APP_PSEN1_PSEN2_common_drugTargets.promoter.down.DEGs <- merge(APP_PSEN1_common_drugTargets.promoter.down.DEGs, drugTargets_dATACpeaks.promoter.N141I.DEG.down, by = "Drug_IDs")


#Write .csv files
write.csv(APP_PSEN1_common_drugTargets.promoter.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Down.DEGs/230803_Common_Drug_Targets_Promoter.Down_APP_PSEN1.csv", sep="\t", quote=F, row.names=T)
write.csv(PSEN1_PSEN2_common_drugTargets.promoter.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Down.DEGs/230803_Common_Drug_Targets_Promoter.Down_PSEN1_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN2_common_drugTargets.promoter.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Down.DEGs/230803_Common_Drug_Targets_Promoter.Down_APP_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN1_PSEN2_common_drugTargets.promoter.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Down.DEGs/230803_Common_Drug_Targets_Promoter.Down_allMutations.csv", sep="\t", quote=F, row.names=T)

#################################################################
#35. Create Venn Diagrams for drug targets based on integrated
#DEGs and dATAC peaks - promoter down regions
#################################################################

#Create drug target list for APP-V717I 
drugTargets_V717I.promoter.down_list <- list(drugTargets_dATACpeaks.promoter.V717I.DEG.down = as.character(unique(drugTargets_dATACpeaks.promoter.V717I.DEG.down$Drug_IDs)))

#Create drug target list for PSEN1-A79V 
drugTargets_A79V.promoter.down_list <- list(drugTargets_dATACpeaks.promoter.A79V.DEG.down = as.character(unique(drugTargets_dATACpeaks.promoter.A79V.DEG.down$Drug_IDs)))

#Create drug target list for PSEN2-N141I 
drugTargets_N141I.promoter.down_list <- list(drugTargets_dATACpeaks.promoter.N141I.DEG.down = as.character(unique(drugTargets_dATACpeaks.promoter.N141I.DEG.down$Drug_IDs)))

#Plot the nVenn Diagram 
myNV_dATAC_DEGs_drugTargets_promoter.down <- plotVenn(list(drugTargets_V717I.promoter.down_list, drugTargets_A79V.promoter.down_list, drugTargets_N141I.promoter.down_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)

#Show the .svg file
showSVG(myNV_dATAC_DEGs_drugTargets_promoter.down, opacity=0.3,outFile = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_promoter.down.DEGs_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))


######################################################################
#36. Create Chord Diagram for subset of drug targets across all mutations
#Promoter Down Regions
######################################################################

#Get unique gene symbols across common drug targets of all mutations
#to make input matrix with 0's and 1's
chosenDrugTargets.Promoter.Down <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Down.DEGs/230803_Chosen_Drug_Targets_Promoter.Down_allMutations.csv")

#Filter by drug and get unique symbols
uniqueGenes.DB01254.Promoter.Down1 <- chosenDrugTargets.Promoter.Down %>%
  dplyr::filter(Drug_IDs == "DB01254") 
uniqueGenes.DB01254.Promoter.Down.1a <- data.frame(unique(uniqueGenes.DB01254.Promoter.Down1$SYMBOL.x))
uniqueGenes.DB01254.Promoter.Down.1b <- data.frame(unique(uniqueGenes.DB01254.Promoter.Down1$SYMBOL.y))
uniqueGenes.DB01254.Promoter.Down.1c <- data.frame(unique(uniqueGenes.DB01254.Promoter.Down1$SYMBOL))

uniqueGenes.DB01268.Promoter.Down2 <- chosenDrugTargets.Promoter.Down %>%
  dplyr::filter(Drug_IDs == "DB01268") 
uniqueGenes.DB01268.Promoter.Down.2a <- data.frame(unique(uniqueGenes.DB01268.Promoter.Down2$SYMBOL.x))
uniqueGenes.DB01268.Promoter.Down.2b <- data.frame(unique(uniqueGenes.DB01268.Promoter.Down2$SYMBOL.y))
uniqueGenes.DB01268.Promoter.Down.2c <- data.frame(unique(uniqueGenes.DB01268.Promoter.Down2$SYMBOL))

uniqueGenes.DB06616.Promoter.Down3 <- chosenDrugTargets.Promoter.Down %>%
  dplyr::filter(Drug_IDs == "DB06616") 
uniqueGenes.DB06616.Promoter.Down.3a <- data.frame(unique(uniqueGenes.DB06616.Promoter.Down3$SYMBOL.x))
uniqueGenes.DB06616.Promoter.Down.3b <- data.frame(unique(uniqueGenes.DB06616.Promoter.Down3$SYMBOL.y))
uniqueGenes.DB06616.Promoter.Down.3c <- data.frame(unique(uniqueGenes.DB06616.Promoter.Down3$SYMBOL))

uniqueGenes.DB08865.Promoter.Down4 <- chosenDrugTargets.Promoter.Down %>%
  dplyr::filter(Drug_IDs == "DB08865") 
uniqueGenes.DB08865.Promoter.Down.4a <- data.frame(unique(uniqueGenes.DB08865.Promoter.Down4$SYMBOL.x))
uniqueGenes.DB08865.Promoter.Down.4b <- data.frame(unique(uniqueGenes.DB08865.Promoter.Down4$SYMBOL.y))
uniqueGenes.DB08865.Promoter.Down.4c <- data.frame(unique(uniqueGenes.DB08865.Promoter.Down4$SYMBOL))

uniqueGenes.DB08877.Promoter.Down5 <- chosenDrugTargets.Promoter.Down %>%
  dplyr::filter(Drug_IDs == "DB08877") 
uniqueGenes.DB08877.Promoter.Down.5a <- data.frame(unique(uniqueGenes.DB08877.Promoter.Down5$SYMBOL.x))
uniqueGenes.DB08877.Promoter.Down.5b <- data.frame(unique(uniqueGenes.DB08877.Promoter.Down5$SYMBOL.y))
uniqueGenes.DB08877.Promoter.Down.5c <- data.frame(unique(uniqueGenes.DB08877.Promoter.Down5$SYMBOL))

uniqueGenes.DB09079.Promoter.Down6 <- chosenDrugTargets.Promoter.Down %>%
  dplyr::filter(Drug_IDs == "DB09079") 
uniqueGenes.DB09079.Promoter.Down.6a <- data.frame(unique(uniqueGenes.DB09079.Promoter.Down6$SYMBOL.x))
uniqueGenes.DB09079.Promoter.Down.6b <- data.frame(unique(uniqueGenes.DB09079.Promoter.Down6$SYMBOL.y))
uniqueGenes.DB09079.Promoter.Down.6c <- data.frame(unique(uniqueGenes.DB09079.Promoter.Down6$SYMBOL))

#Read .csv file with drug values across PSEN1/PSEN2/APP
chordDiagram_Values_Promoter.Down <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/Promoter_Down.DEGs/230803_Promoter.Down_chordDiag.Values_allMutations.csv")

#Source: https://urldefense.com/v3/__https://stackoverflow.com/questions/73609132/how-to-create-a-chord-diagram-in-r__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqcD_PAlf$ 
#Set the colors for each drug
#Source: https://urldefense.com/v3/__https://r-charts.com/flow/chord-diagram/__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqe_Um70Q$ 
colors2 <- c(Dasatinib = "#009e73",Sunitinib = "#449AE4", 
             Bosutinib = "#871C9A", Crizotinib = "#2C2E9C",
             Ruxolitinib = "#820D3F", Nintedanib = "#E64A00")

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Promoter.Down_Chord_labeled.svg")

set.seed(123)
chordDiagram(chordDiagram_Values_Promoter.Down, grid.col = colors2)
dev.off()

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Promoter.Down_Chord_unlabeled.svg")

set.seed(123)
chordDiagram(chordDiagram_Values_Promoter.Down, grid.col = colors2, annotationTrack = "grid")
dev.off()


#####################################################################
#37. Get the common drug targets by ID for PERGRINE enhancer up regions
#####################################################################

#Get common drug targets by drugID
APP_PSEN1_common_drugTargets.enhancer.PEREGRINE.up.DEGs <- merge(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.up, drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.up, by = "Drug_IDs")
PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.up.DEGs <- merge(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.up, drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.up, by = "Drug_IDs")
APP_PSEN2_common_drugTargets.enhancer.PEREGRINE.up.DEGs <- merge(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.up, drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.up, by = "Drug_IDs")
APP_PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.up.DEGs <- merge(APP_PSEN1_common_drugTargets.enhancer.PEREGRINE.up.DEGs, drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.up, by = "Drug_IDs")


#Write .csv files
write.csv(APP_PSEN1_common_drugTargets.enhancer.PEREGRINE.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Up.DEGs/230803_Common_Drug_Targets_enhancer.PEREGRINE.Up_APP_PSEN1.csv", sep="\t", quote=F, row.names=T)
write.csv(PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Up.DEGs/230803_Common_Drug_Targets_enhancer.PEREGRINE.Up_PSEN1_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN2_common_drugTargets.enhancer.PEREGRINE.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Up.DEGs/230803_Common_Drug_Targets_enhancer.PEREGRINE.Up_APP_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.table(APP_PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.up.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Up.DEGs/230804_Common_Drug_Targets_enhancer.PEREGRINE.Up_allMutations.txt", sep="\t", quote=F, row.names=T)


#################################################################
#38. Create Venn Diagram for drug targets based on integrated
#DEGs and dATAC peaks - PEREGRINE enhancer up regions
#################################################################

#Create drug target list for APP-V717I 
drugTargets_V717I.enhancer.PEREGRINE.up_list <- list(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.up = as.character(unique(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.up$Drug_IDs)))

#Create drug target list for PSEN1-A79V 
drugTargets_A79V.enhancer.PEREGRINE.up_list <- list(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.up = as.character(unique(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.up$Drug_IDs)))

#Create drug target list for PSEN2-N141I 
drugTargets_N141I.enhancer.PEREGRINE.up_list <- list(drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.up = as.character(unique(drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.up$Drug_IDs)))

#Plot the nVenn Diagram 
myNV_dATAC_DEGs_drugTargets_enhancer.PEREGRINE.up <- plotVenn(list(drugTargets_V717I.enhancer.PEREGRINE.up_list, drugTargets_A79V.enhancer.PEREGRINE.up_list, drugTargets_N141I.enhancer.PEREGRINE.up_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)

#Show the .svg file
showSVG(myNV_dATAC_DEGs_drugTargets_enhancer.PEREGRINE.up, opacity=0.3,outFile = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_enhancer.PEREGRINE.up.DEGs_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))


######################################################################
#39. Create Chord Diagram for subset of drug targets across all mutations
#PEERGRINE Enhancer Up Regions
######################################################################

#Get unique gene symbols across common drug targets of all mutations
#to make input matrix with 0's and 1's
chosenDrugTargets.PEREGRINE.Enhancer.Up <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Up.DEGs/230804_Chosen_Drug_Targets_PEREGRINE.Enhancer.Up_allMutations.csv")

#Filter by drug and get unique symbols
uniqueGenes.DB00898.PEREGRINE.Enhancer.Up1 <- chosenDrugTargets.PEREGRINE.Enhancer.Up %>%
  dplyr::filter(Drug_IDs == "DB00898") 
uniqueGenes.DB00898.PEREGRINE.Enhancer.Up.1a <- data.frame(unique(uniqueGenes.DB00898.PEREGRINE.Enhancer.Up1$SYMBOL.x))
uniqueGenes.DB00898.PEREGRINE.Enhancer.Up.1b <- data.frame(unique(uniqueGenes.DB00898.PEREGRINE.Enhancer.Up1$SYMBOL.y))
uniqueGenes.DB00898.PEREGRINE.Enhancer.Up.1c <- data.frame(unique(uniqueGenes.DB00898.PEREGRINE.Enhancer.Up1$SYMBOL))

uniqueGenes.DB01268.PEREGRINE.Enhancer.Up2 <- chosenDrugTargets.PEREGRINE.Enhancer.Up %>%
  dplyr::filter(Drug_IDs == "DB01268") 
uniqueGenes.DB01268.PEREGRINE.Enhancer.Up.2a <- data.frame(unique(uniqueGenes.DB01268.PEREGRINE.Enhancer.Up2$SYMBOL.x))
uniqueGenes.DB01268.PEREGRINE.Enhancer.Up.2b <- data.frame(unique(uniqueGenes.DB01268.PEREGRINE.Enhancer.Up2$SYMBOL.y))
uniqueGenes.DB01268.PEREGRINE.Enhancer.Up.2c <- data.frame(unique(uniqueGenes.DB01268.PEREGRINE.Enhancer.Up2$SYMBOL))

uniqueGenes.DB06616.PEREGRINE.Enhancer.Up3 <- chosenDrugTargets.PEREGRINE.Enhancer.Up %>%
  dplyr::filter(Drug_IDs == "DB06616") 
uniqueGenes.DB06616.PEREGRINE.Enhancer.Up.3a <- data.frame(unique(uniqueGenes.DB06616.PEREGRINE.Enhancer.Up3$SYMBOL.x))
uniqueGenes.DB06616.PEREGRINE.Enhancer.Up.3b <- data.frame(unique(uniqueGenes.DB06616.PEREGRINE.Enhancer.Up3$SYMBOL.y))
uniqueGenes.DB06616.PEREGRINE.Enhancer.Up.3c <- data.frame(unique(uniqueGenes.DB06616.PEREGRINE.Enhancer.Up3$SYMBOL))


uniqueGenes.DB08865.PEREGRINE.Enhancer.Up4 <- chosenDrugTargets.PEREGRINE.Enhancer.Up %>%
  dplyr::filter(Drug_IDs == "DB08865") 
uniqueGenes.DB08865.PEREGRINE.Enhancer.Up.4a <- data.frame(unique(uniqueGenes.DB08865.PEREGRINE.Enhancer.Up4$SYMBOL.x))
uniqueGenes.DB08865.PEREGRINE.Enhancer.Up.4b <- data.frame(unique(uniqueGenes.DB08865.PEREGRINE.Enhancer.Up4$SYMBOL.y))
uniqueGenes.DB08865.PEREGRINE.Enhancer.Up.4c <- data.frame(unique(uniqueGenes.DB08865.PEREGRINE.Enhancer.Up4$SYMBOL))


uniqueGenes.DB08877.PEREGRINE.Enhancer.Up5 <- chosenDrugTargets.PEREGRINE.Enhancer.Up %>%
  dplyr::filter(Drug_IDs == "DB08877") 
uniqueGenes.DB08877.PEREGRINE.Enhancer.Up.5a <- data.frame(unique(uniqueGenes.DB08877.PEREGRINE.Enhancer.Up5$SYMBOL.x))
uniqueGenes.DB08877.PEREGRINE.Enhancer.Up.5b <- data.frame(unique(uniqueGenes.DB08877.PEREGRINE.Enhancer.Up5$SYMBOL.y))
uniqueGenes.DB08877.PEREGRINE.Enhancer.Up.5c <- data.frame(unique(uniqueGenes.DB08877.PEREGRINE.Enhancer.Up5$SYMBOL))

uniqueGenes.DB09079.PEREGRINE.Enhancer.Up6 <- chosenDrugTargets.PEREGRINE.Enhancer.Up %>%
  dplyr::filter(Drug_IDs == "DB09079") 
uniqueGenes.DB09079.PEREGRINE.Enhancer.Up.6a <- data.frame(unique(uniqueGenes.DB09079.PEREGRINE.Enhancer.Up6$SYMBOL.x))
uniqueGenes.DB09079.PEREGRINE.Enhancer.Up.6b <- data.frame(unique(uniqueGenes.DB09079.PEREGRINE.Enhancer.Up6$SYMBOL.y))
uniqueGenes.DB09079.PEREGRINE.Enhancer.Up.6c <- data.frame(unique(uniqueGenes.DB09079.PEREGRINE.Enhancer.Up6$SYMBOL))

uniqueGenes.DB00398.PEREGRINE.Enhancer.Up7 <- chosenDrugTargets.PEREGRINE.Enhancer.Up %>%
  dplyr::filter(Drug_IDs == "DB00398") 
uniqueGenes.DB00398.PEREGRINE.Enhancer.Up.7a <- data.frame(unique(uniqueGenes.DB00398.PEREGRINE.Enhancer.Up7$SYMBOL.x))
uniqueGenes.DB00398.PEREGRINE.Enhancer.Up.7b <- data.frame(unique(uniqueGenes.DB00398.PEREGRINE.Enhancer.Up7$SYMBOL.y))
uniqueGenes.DB00398.PEREGRINE.Enhancer.Up.7c <- data.frame(unique(uniqueGenes.DB00398.PEREGRINE.Enhancer.Up7$SYMBOL))

#Read .csv file with drug values across PSEN1/PSEN2/APP
chordDiagram_Values_PEREGRINE.Enhancer.Up <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Up.DEGs/230804_PEREGRINE.Enhancer.Up_chordDiag.Values_allMutations.csv")

#Source: https://urldefense.com/v3/__https://stackoverflow.com/questions/73609132/how-to-create-a-chord-diagram-in-r__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqcD_PAlf$ 
#Set the colors for each drug
#Source: https://urldefense.com/v3/__https://r-charts.com/flow/chord-diagram/__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqe_Um70Q$ 
colors3 <- c(Ethanol = "#FFC107", Sunitinib = "#449AE4", 
             Bosutinib = "#871C9A", Sorafenib = "#CC79A7", 
             Ruxolitinib = "#820D3F", Crizotinib = "#2C2E9C", 
             Nintedanib = "#E64A00"
  )

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_PEREGRINE.Enhancer.Up_Chord_labeled.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_PEREGRINE.Enhancer.Up, grid.col = colors3)
dev.off()

#Output to an .svg file and create the unlabeled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_PEREGRINE.Enhancer.Up_Chord_unlabeled.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_PEREGRINE.Enhancer.Up, grid.col = colors3, annotationTrack = "grid")
dev.off()


#####################################################################
#40. Get the common drug targets by ID for PEREGRINE enhancer down regions
#####################################################################

#Get common drug targets by drugID
APP_PSEN1_common_drugTargets.enhancer.PEREGRINE.down.DEGs <- merge(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.down, drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.down, by = "Drug_IDs")
PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.down.DEGs <- merge(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.down, drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.down, by = "Drug_IDs")
APP_PSEN2_common_drugTargets.enhancer.PEREGRINE.down.DEGs <- merge(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.down, drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.down, by = "Drug_IDs")
APP_PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.down.DEGs <- merge(APP_PSEN1_common_drugTargets.enhancer.PEREGRINE.down.DEGs, drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.down, by = "Drug_IDs")


#Write .csv files
write.csv(APP_PSEN1_common_drugTargets.enhancer.PEREGRINE.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Down.DEGs/230803_Common_Drug_Targets_enhancer.PEREGRINE.Down_APP_PSEN1.csv", sep="\t", quote=F, row.names=T)
write.csv(PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Down.DEGs/230803_Common_Drug_Targets_enhancer.PEREGRINE.Down_PSEN1_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN2_common_drugTargets.enhancer.PEREGRINE.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Down.DEGs/230803_Common_Drug_Targets_enhancer.PEREGRINE.Down_APP_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.table(APP_PSEN1_PSEN2_common_drugTargets.enhancer.PEREGRINE.down.DEGs, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Down.DEGs/230804_Common_Drug_Targets_enhancer.PEREGRINE.Down_allMutations.txt", sep="\t", quote=F, row.names=T)

#################################################################
#41. Create Venn Diagram for drug targets based on integrated
#DEGs and dATAC peaks - PEREGRINE enhancer down regions
#################################################################

#Create drug target list for APP-V717I 
drugTargets_V717I.enhancer.PEREGRINE.down_list <- list(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.down = as.character(unique(drugTargets_dATACpeaks.enhancer.PEREGRINE.V717I.DEG.down$Drug_IDs)))

#Create drug target list for PSEN1-A79V 
drugTargets_A79V.enhancer.PEREGRINE.down_list <- list(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.down = as.character(unique(drugTargets_dATACpeaks.enhancer.PEREGRINE.A79V.DEG.down$Drug_IDs)))

#Create drug target list for PSEN2-N141I 
drugTargets_N141I.enhancer.PEREGRINE.down_list <- list(drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.down = as.character(unique(drugTargets_dATACpeaks.enhancer.PEREGRINE.N141I.DEG.down$Drug_IDs)))

#Plot the nVenn Diagram 
myNV_dATAC_DEGs_drugTargets_enhancer.PEREGRINE.down <- plotVenn(list(drugTargets_V717I.enhancer.PEREGRINE.down_list, drugTargets_A79V.enhancer.PEREGRINE.up_list, drugTargets_N141I.enhancer.PEREGRINE.down_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)

#Show the .svg file
showSVG(myNV_dATAC_DEGs_drugTargets_enhancer.PEREGRINE.down, opacity=0.3,outFile = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_enhancer.PEREGRINE.down.DEGs_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))


######################################################################
#42. Create Chord Diagram for subset of drug targets across all mutations
#PEERGRINE Enhancer Down Regions
######################################################################

#Get unique gene symbols across common drug targets of all mutations
#to make input matrix with 0's and 1's
chosenDrugTargets.PEREGRINE.Enhancer.Down <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Down.DEGs/230804_Chosen_Drug_Targets_PEREGRINE.Enhancer.Down_allMutations.csv")

#Filter by drug and get unique symbols
uniqueGenes.DB00996.PEREGRINE.Enhancer.Down1 <- chosenDrugTargets.PEREGRINE.Enhancer.Down %>%
  dplyr::filter(Drug_IDs == "DB00996") 
uniqueGenes.DB00996.PEREGRINE.Enhancer.Down.1a <- data.frame(unique(uniqueGenes.DB00996.PEREGRINE.Enhancer.Down1$SYMBOL.x))
uniqueGenes.DB00996.PEREGRINE.Enhancer.Down.1b <- data.frame(unique(uniqueGenes.DB00996.PEREGRINE.Enhancer.Down1$SYMBOL.y))
uniqueGenes.DB00996.PEREGRINE.Enhancer.Down.1c <- data.frame(unique(uniqueGenes.DB00996.PEREGRINE.Enhancer.Down1$SYMBOL))

uniqueGenes.DB01254.PEREGRINE.Enhancer.Down2 <- chosenDrugTargets.PEREGRINE.Enhancer.Down %>%
  dplyr::filter(Drug_IDs == "DB01254") 
uniqueGenes.DB01254.PEREGRINE.Enhancer.Down.2a <- data.frame(unique(uniqueGenes.DB01254.PEREGRINE.Enhancer.Down2$SYMBOL.x))
uniqueGenes.DB01254.PEREGRINE.Enhancer.Down.2b <- data.frame(unique(uniqueGenes.DB01254.PEREGRINE.Enhancer.Down2$SYMBOL.y))
uniqueGenes.DB01254.PEREGRINE.Enhancer.Down.2c <- data.frame(unique(uniqueGenes.DB01254.PEREGRINE.Enhancer.Down2$SYMBOL))

uniqueGenes.DB01268.PEREGRINE.Enhancer.Down3 <- chosenDrugTargets.PEREGRINE.Enhancer.Down %>%
  dplyr::filter(Drug_IDs == "DB01268") 
uniqueGenes.DB01268.PEREGRINE.Enhancer.Down.3a <- data.frame(unique(uniqueGenes.DB01268.PEREGRINE.Enhancer.Down3$SYMBOL.x))
uniqueGenes.DB01268.PEREGRINE.Enhancer.Down.3b <- data.frame(unique(uniqueGenes.DB01268.PEREGRINE.Enhancer.Down3$SYMBOL.y))
uniqueGenes.DB01268.PEREGRINE.Enhancer.Down.3c <- data.frame(unique(uniqueGenes.DB01268.PEREGRINE.Enhancer.Down3$SYMBOL))


uniqueGenes.DB04868.PEREGRINE.Enhancer.Down4 <- chosenDrugTargets.PEREGRINE.Enhancer.Down %>%
  dplyr::filter(Drug_IDs == "DB04868") 
uniqueGenes.DB04868.PEREGRINE.Enhancer.Down.4a <- data.frame(unique(uniqueGenes.DB04868.PEREGRINE.Enhancer.Down4$SYMBOL.x))
uniqueGenes.DB04868.PEREGRINE.Enhancer.Down.4b <- data.frame(unique(uniqueGenes.DB04868.PEREGRINE.Enhancer.Down4$SYMBOL.y))
uniqueGenes.DB04868.PEREGRINE.Enhancer.Down.4c <- data.frame(unique(uniqueGenes.DB04868.PEREGRINE.Enhancer.Down4$SYMBOL))


uniqueGenes.DB06616.PEREGRINE.Enhancer.Down5 <- chosenDrugTargets.PEREGRINE.Enhancer.Down %>%
  dplyr::filter(Drug_IDs == "DB06616") 
uniqueGenes.DB06616.PEREGRINE.Enhancer.Down.5a <- data.frame(unique(uniqueGenes.DB06616.PEREGRINE.Enhancer.Down5$SYMBOL.x))
uniqueGenes.DB06616.PEREGRINE.Enhancer.Down.5b <- data.frame(unique(uniqueGenes.DB06616.PEREGRINE.Enhancer.Down5$SYMBOL.y))
uniqueGenes.DB06616.PEREGRINE.Enhancer.Down.5c <- data.frame(unique(uniqueGenes.DB06616.PEREGRINE.Enhancer.Down5$SYMBOL))

uniqueGenes.DB08877.PEREGRINE.Enhancer.Down6 <- chosenDrugTargets.PEREGRINE.Enhancer.Down %>%
  dplyr::filter(Drug_IDs == "DB08877") 
uniqueGenes.DB08877.PEREGRINE.Enhancer.Down.6a <- data.frame(unique(uniqueGenes.DB08877.PEREGRINE.Enhancer.Down6$SYMBOL.x))
uniqueGenes.DB08877.PEREGRINE.Enhancer.Down.6b <- data.frame(unique(uniqueGenes.DB08877.PEREGRINE.Enhancer.Down6$SYMBOL.y))
uniqueGenes.DB08877.PEREGRINE.Enhancer.Down.6c <- data.frame(unique(uniqueGenes.DB08877.PEREGRINE.Enhancer.Down6$SYMBOL))

uniqueGenes.DB01183.PEREGRINE.Enhancer.Down7 <- chosenDrugTargets.PEREGRINE.Enhancer.Down %>%
  dplyr::filter(Drug_IDs == "DB01183") 
uniqueGenes.DB01183.PEREGRINE.Enhancer.Down.7a <- data.frame(unique(uniqueGenes.DB01183.PEREGRINE.Enhancer.Down7$SYMBOL.x))
uniqueGenes.DB01183.PEREGRINE.Enhancer.Down.7b <- data.frame(unique(uniqueGenes.DB01183.PEREGRINE.Enhancer.Down7$SYMBOL.y))
uniqueGenes.DB01183.PEREGRINE.Enhancer.Down.7c <- data.frame(unique(uniqueGenes.DB01183.PEREGRINE.Enhancer.Down7$SYMBOL))

#Read .csv file with drug input matrix of 0's and 1's
chordDiagram_Values_PEREGRINE.Enhancer.Down <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/PEREGRINE.Enhancer_Down.DEGs/230804_PEREGRINE.Enhancer.Down_chordDiag.Values_allMutations.csv")

#Source: https://urldefense.com/v3/__https://stackoverflow.com/questions/73609132/how-to-create-a-chord-diagram-in-r__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqcD_PAlf$ 
#Set the colors for each drug
#Source: https://urldefense.com/v3/__https://r-charts.com/flow/chord-diagram/__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqe_Um70Q$ 
colors4 <- c(Gabapentin = "#e69f00", Dasatinib = "#009e73",
             Sunitinib = "#449AE4", Nilotinib = "#14C7BA",
             Bosutinib = "#871C9A", Ruxolitinib = "#820D3F",
             Naloxone = "#E884AF" 
)

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_PEREGRINE.Enhancer.Down_Chord_labeled.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_PEREGRINE.Enhancer.Down, grid.col = colors4)
dev.off()

#Output to an .svg file and create the unlabeled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_PEREGRINE.Enhancer.Down_Chord_unlabeled.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_PEREGRINE.Enhancer.Down, grid.col = colors4, annotationTrack = "grid")
dev.off()

##=====================================================================================================================================================================================================================================================================================

########################################################################
#43. Get the common drug targets from all phase trials by CERNO pathways 
#ranked by intePareto for all mutations
########################################################################

APP_PSEN1_common_drugTargets.CERNO.CADROclasses <- merge(CERNO_CADROclasses.drugTargets.V717I, CERNO_CADROclasses.drugTargets.A79V, by = "Title")
PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses <- merge(CERNO_CADROclasses.drugTargets.A79V, CERNO_CADROclasses.drugTargets.N141I, by = "Title")
APP_PSEN2_common_drugTargets.CERNO.CADROclasses <- merge(CERNO_CADROclasses.drugTargets.V717I, CERNO_CADROclasses.drugTargets.N141I, by = "Title")
APP_PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses <- merge(APP_PSEN1_common_drugTargets.CERNO.CADROclasses, CERNO_CADROclasses.drugTargets.N141I, by = "Title")
 
# #Write.csv files
write.csv(APP_PSEN1_common_drugTargets.CERNO.CADROclasses, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/APP_PSEN1_common_drugTargets.CERNO.CADROclasses.csv", sep="\t", quote=F, row.names=T)
write.csv(PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses.csv", sep="\t", quote=F, row.names=T) write.csv(APP_PSEN2_common_drugTargets.CERNO.CADROclasses, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/APP_PSEN2_common_drugTargets.CERNO.CADROclasses.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/APP_PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses.csv", sep="\t", quote=F, row.names=T)
write.table(APP_PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses, file="/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/APP_PSEN1_PSEN2_common_drugTargets.CERNO.CADROclasses.txt", sep="\t", quote=F, row.names=T)



#################################################################
#44. Create Venn Diagram for CADRO drug agents for Phase Trial 1
#based on CERNO pathways ranked by intePareto
#################################################################

#Subset to only Phase 1 trial drug agents for APP-V717I
CERNO_CADROclasses.V717I.Phase1 <- CERNO_CADROclasses.drugTargets.V717I %>%
  dplyr::filter(Phase_Trial == "Phase 1")

#Subset to only Phase 1 trial drug agents for PSEN1-A79V
CERNO_CADROclasses.A79V.Phase1 <- CERNO_CADROclasses.drugTargets.A79V %>%
  dplyr::filter(Phase_Trial == "Phase 1")

#Subset to only Phase 1 trial drug agents for PSEN2-N141I
CERNO_CADROclasses.N141I.Phase1 <- CERNO_CADROclasses.drugTargets.N141I %>%
  dplyr::filter(Phase_Trial == "Phase 1")

#Create CADRO drug agent list for APP-V717I 
CADRO_drugTargets.Phase1_V717I_list <- list(CERNO_CADROclasses.V717I.Phase1  = as.character(unique(CERNO_CADROclasses.V717I.Phase1$Agent)))

#Create CADRO drug agent list for PSEN1-A79V 
CADRO_drugTargets.Phase1_A79V_list <- list(CERNO_CADROclasses.A79V.Phase1  = as.character(unique(CERNO_CADROclasses.A79V.Phase1$Agent)))

#Create CADRO drug agent list for PSEN2-N141I 
CADRO_drugTargets.Phase1_N141I_list <- list(CERNO_CADROclasses.N141I.Phase1  = as.character(unique(CERNO_CADROclasses.N141I.Phase1$Agent)))

#Plot the nVenn Diagram 
myNV_CADRO_drugTargets.Phase1<- plotVenn(list(CADRO_drugTargets.Phase1_V717I_list, CADRO_drugTargets.Phase1_A79V_list, CADRO_drugTargets.Phase1_N141I_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)

#Show the .svg file
showSVG(myNV_CADRO_drugTargets.Phase1, opacity=0.3,outFile = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_CERNO.CADRO_drugTargets.Phase1_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))


######################################################################
#45. Create Chord Diagram for subset of Phase 1 drug agents with defined
#CADRO mechanism classes and pathways across all mutations
######################################################################

#Get unique drug agents across common enriched pathways of all mutations
#to make input matrix with 0's and 1's
chosenDrugPathways.CERNO.CADROclasses <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/230804_chosen_drugTargets.CERNO.CADROclasses_allMutations.csv")

#Filter by pathway and drug trial phase and then get unique drug agents
#to create matrix of 0's and 1's for input matrix

#Biological Adhesion
uniqueDrugPathway.biologicalAdhesion.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp biological adhesion") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.biologicalAdhesion.Phase1 <- data.frame(unique(uniqueDrugPathway.biologicalAdhesion.Phase1$Agent))

#Cell Fate Specification
uniqueDrugPathway.cellFate.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp cell fate commitment") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.cellFate.Phase1 <- data.frame(unique(uniqueDrugPathway.cellFate.Phase1$Agent))

#CNS Neuron Differentiation
uniqueDrugPathway.CNS_neuronDiff.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp central nervous system neuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.CNS_neuronDiff.Phase1 <- data.frame(unique(uniqueDrugPathway.CNS_neuronDiff.Phase1$Agent))

#Cerebral Cortex GABAergic Interneuron Differentiation
uniqueDrugPathway.cortex_GABAergic.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp cerebral cortex gabaergic interneuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.cortex_GABAergic.Phase1 <- data.frame(unique(uniqueDrugPathway.cortex_GABAergic.Phase1$Agent))

#G protein coupled receptor signaling
uniqueDrugPathway.Gprotein_signaling.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp g protein coupled receptor signaling pathway") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.Gprotein_signaling.Phase1 <- data.frame(unique(uniqueDrugPathway.Gprotein_signaling.Phase1$Agent))

#Nervous system process
uniqueDrugPathway.nervousSystem.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp nervous system process") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.nervousSystem.Phase1 <- data.frame(unique(uniqueDrugPathway.nervousSystem.Phase1$Agent))

#Neuron Fate Committment
uniqueDrugPathway.neuronFateCommit.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp neuron fate commitment") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.neuronFateCommit.Phase1 <- data.frame(unique(uniqueDrugPathway.neuronFateCommit.Phase1$Agent))

#Taxis
uniqueDrugPathway.taxis.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp taxis") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.taxis.Phase1 <- data.frame(unique(uniqueDrugPathway.taxis.Phase1$Agent))

#Neuron Fate Specification
uniqueDrugPathway.neuronFateSpec.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp neuron fate specification") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.neuronFateSpec.Phase1 <- data.frame(unique(uniqueDrugPathway.neuronFateSpec.Phase1$Agent))

#Regulation of neuron differentiation
uniqueDrugPathway.reg_neuronDiff.Phase1 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp regulation of neuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 1")
uniqueDrugAgent.reg_neuronDiff.Phase1 <- data.frame(unique(uniqueDrugPathway.reg_neuronDiff.Phase1$Agent))

#Read .csv file with drug input matrix of 0's and 1's
chordDiagram_Values_drugTargets.Phase1 <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/230805_chordDiagValues_allMutations_drugTargets_Phase1.csv")

#Source: https://urldefense.com/v3/__https://stackoverflow.com/questions/73609132/how-to-create-a-chord-diagram-in-r__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqcD_PAlf$ 
#Set the colors for each drug
#Source: https://urldefense.com/v3/__https://r-charts.com/flow/chord-diagram/__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqe_Um70Q$ 
colors5.Phase1 <- c(BDPP = "#004949",
                    BEY2153 = "#78B389", 
                    Contraloid.acetate = "#008800",
                    Allopregnanolone = "#D35FB7",
                    NNI.362 = "#ffb6db",
                    BMS.984923 = "#490092",
                    COR588 = "#b66dff",
                    IGC.AD1 = "#006ddb",
                    MK1942.donepezil = "#6db6ff",
                    Donepezil = "#b6dbff",
                    Edicotinib = "#882255",
                    Salsalate = "#FFC107",
                    VT301 = "#FFC20A",
                    XPro1595 = "#F1BB83",
                    Emtricitabine = "#E1BE6A" 
)

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Phase1_Chord_labeled_v2.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_drugTargets.Phase1, grid.col = colors5.Phase1)
dev.off()

#Output to an .svg file and create the unlabeled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Phase1_Chord_unlabeled_v2.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_drugTargets.Phase1, grid.col = colors5.Phase1, annotationTrack = "grid")
dev.off()


#################################################################
#46. Create Venn Diagram for CADRO drug agents for Phase Trial 2
#based on CERNO pathways ranked by intePareto
#################################################################

#Subset to only Phase 2 trial drug agents for APP-V717I
CERNO_CADROclasses.V717I.Phase2 <- CERNO_CADROclasses.drugTargets.V717I %>%
  dplyr::filter(Phase_Trial == "Phase 2")

#Subset to only Phase 2 trial drug agents for PSEN1-A79V
CERNO_CADROclasses.A79V.Phase2 <- CERNO_CADROclasses.drugTargets.A79V %>%
  dplyr::filter(Phase_Trial == "Phase 2")

#Subset to only Phase 2 trial drug agents for PSEN2-N141I
CERNO_CADROclasses.N141I.Phase2 <- CERNO_CADROclasses.drugTargets.N141I %>%
  dplyr::filter(Phase_Trial == "Phase 2")

#Create CADRO drug agent list for APP-V717I 
CADRO_drugTargets.Phase2_V717I_list <- list(CERNO_CADROclasses.V717I.Phase2  = as.character(unique(CERNO_CADROclasses.V717I.Phase2$Agent)))

#Create CADRO drug agent list for PSEN1-A79V 
CADRO_drugTargets.Phase2_A79V_list <- list(CERNO_CADROclasses.A79V.Phase2  = as.character(unique(CERNO_CADROclasses.A79V.Phase2$Agent)))

#Create CADRO drug agent list for PSEN2-N141I 
CADRO_drugTargets.Phase2_N141I_list <- list(CERNO_CADROclasses.N141I.Phase2  = as.character(unique(CERNO_CADROclasses.N141I.Phase2$Agent)))

#Plot the nVenn Diagram 
myNV_CADRO_drugTargets.Phase2<- plotVenn(list(CADRO_drugTargets.Phase2_V717I_list, CADRO_drugTargets.Phase2_A79V_list, CADRO_drugTargets.Phase2_N141I_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)

#Show the .svg file
showSVG(myNV_CADRO_drugTargets.Phase2, opacity=0.3,outFile = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_CERNO.CADRO_drugTargets.Phase2_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))


######################################################################
#47. Create Chord Diagram for subset of Phase 2 drug agents with defined
#CADRO mechanism classes and pathways across all mutations
######################################################################

#Filter by pathway and drug trial phase and then get unique drug agents
#to create matrix of 0's and 1's for input matrix

#Biological Adhesion
uniqueDrugPathway.biologicalAdhesion.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp biological adhesion") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.biologicalAdhesion.Phase2 <- data.frame(unique(uniqueDrugPathway.biologicalAdhesion.Phase2$Agent))

#Cell Fate Specification
uniqueDrugPathway.cellFate.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp cell fate commitment") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.cellFate.Phase2 <- data.frame(unique(uniqueDrugPathway.cellFate.Phase2$Agent))

#CNS Neuron Differentiation
uniqueDrugPathway.CNS_neuronDiff.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp central nervous system neuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.CNS_neuronDiff.Phase2 <- data.frame(unique(uniqueDrugPathway.CNS_neuronDiff.Phase2$Agent))

#Cerebral Cortex GABAergic Interneuron Differentiation
uniqueDrugPathway.cortex_GABAergic.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp cerebral cortex gabaergic interneuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.cortex_GABAergic.Phase2 <- data.frame(unique(uniqueDrugPathway.cortex_GABAergic.Phase2$Agent))

#G protein coupled receptor signaling
uniqueDrugPathway.Gprotein_signaling.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp g protein coupled receptor signaling pathway") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.Gprotein_signaling.Phase2 <- data.frame(unique(uniqueDrugPathway.Gprotein_signaling.Phase2$Agent))

#Nervous system process
uniqueDrugPathway.nervousSystem.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp nervous system process") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.nervousSystem.Phase2 <- data.frame(unique(uniqueDrugPathway.nervousSystem.Phase2$Agent))

#Neuron Fate Committment
uniqueDrugPathway.neuronFateCommit.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp neuron fate commitment") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.neuronFateCommit.Phase2 <- data.frame(unique(uniqueDrugPathway.neuronFateCommit.Phase2$Agent))

#Taxis
uniqueDrugPathway.taxis.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp taxis") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.taxis.Phase2 <- data.frame(unique(uniqueDrugPathway.taxis.Phase2$Agent))

#Neuron Fate Specification
uniqueDrugPathway.neuronFateSpec.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp neuron fate specification") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.neuronFateSpec.Phase2 <- data.frame(unique(uniqueDrugPathway.neuronFateSpec.Phase2$Agent))

#Regulation of neuron differentiation
uniqueDrugPathway.reg_neuronDiff.Phase2 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp regulation of neuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 2")
uniqueDrugAgent.reg_neuronDiff.Phase2 <- data.frame(unique(uniqueDrugPathway.reg_neuronDiff.Phase2$Agent))

#Read .csv file with drug input matrix 
chordDiagram_Values_drugTargets.Phase2 <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/230805_chordDiagValues_allMutations_drugTargets_Phase2.csv")

#Source: https://urldefense.com/v3/__https://stackoverflow.com/questions/73609132/how-to-create-a-chord-diagram-in-r__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqcD_PAlf$ 
#Set the colors for each drug
#Source: https://urldefense.com/v3/__https://r-charts.com/flow/chord-diagram/__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqe_Um70Q$ 
colors6.Phase2 <- c(Posiphen = "#00BB00",
                    Rapamycin = "#117733", 
                    Grapeseed.extract = "#14c7ba",
                    Sovateltide = "#E884AF",
                    Allopregnanolone = "#D35FB7",
                    ExPlas = "#cc79a7",
                    Simufilam = "#B32357",
                    Bryostatin = "#820D3F",
                    Neflamapimod = "#5D3A9B",
                    Elayta = "#4B0092",
                    Edonerpic = "#aa4499",
                    MW150 = "#871c9a",
                    Levetiracetam = "#CE2E68",
                    BPN14770 = "#D41159",
                    Troriluzole = "#905971",
                    Fosgonimeton = "#ff6db6",
                    Dronabinol = "#56b4e9",
                    Bromocriptine = "#009e73",
                    BXCL.501 = "#0072b2",
                    Prazosin = "#009292",
                    AD.35 = "#14c7ba",
                    DAOIB = "#449AE4",
                    Suvorexant = "#3B3EDE",
                    Nicotine.transdermal.patch = "#0C7BDC",
                    Memantine = "#006CD1",
                    THC.free.CBD.oil = "#40B0A6",
                    CST.2032 = "#332288",
                    Baricitinib = "#920000",
                    L.Serine = "#924900",
                    Montelukast = "#db6d00",
                    Pepinemab = "#FD8305",
                    Curcumin.aerobic.yoga = "#E64A00",
                    Canakinumab = "#E66100",
                    TB006 = "#d55e00",
                    AL002 = "#e69f00",
                    Senicapoc = "#994F00"
)

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Phase2_Chord_labeled_v2.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_drugTargets.Phase2, grid.col = colors6.Phase2)
dev.off()

#Output to an .svg file and create the unlabeled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Phase2_Chord_unlabeled_v2.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_drugTargets.Phase2, grid.col = colors6.Phase2, annotationTrack = "grid")
dev.off()



#################################################################
#48. Create Venn Diagram for CADRO drug agents for Phase Trial 3
#based on CERNO pathways ranked by intePareto
#################################################################

#Subset to only Phase 3 trial drug agents for APP-V717I
CERNO_CADROclasses.V717I.Phase3 <- CERNO_CADROclasses.drugTargets.V717I %>%
  dplyr::filter(Phase_Trial == "Phase 3")

#Subset to only Phase 3 trial drug agents for PSEN1-A79V
CERNO_CADROclasses.A79V.Phase3 <- CERNO_CADROclasses.drugTargets.A79V %>%
  dplyr::filter(Phase_Trial == "Phase 3")

#Subset to only Phase 3 trial drug agents for PSEN2-N141I
CERNO_CADROclasses.N141I.Phase3 <- CERNO_CADROclasses.drugTargets.N141I %>%
  dplyr::filter(Phase_Trial == "Phase 3")

#Create CADRO drug agent list for APP-V717I 
CADRO_drugTargets.Phase3_V717I_list <- list(CERNO_CADROclasses.V717I.Phase3  = as.character(unique(CERNO_CADROclasses.V717I.Phase3$Agent)))

#Create CADRO drug agent list for PSEN1-A79V 
CADRO_drugTargets.Phase3_A79V_list <- list(CERNO_CADROclasses.A79V.Phase3  = as.character(unique(CERNO_CADROclasses.A79V.Phase3$Agent)))

#Create CADRO drug agent list for PSEN2-N141I 
CADRO_drugTargets.Phase3_N141I_list <- list(CERNO_CADROclasses.N141I.Phase3  = as.character(unique(CERNO_CADROclasses.N141I.Phase3$Agent)))

#Plot the nVenn Diagram 
myNV_CADRO_drugTargets.Phase3<- plotVenn(list(CADRO_drugTargets.Phase3_V717I_list, CADRO_drugTargets.Phase3_A79V_list, CADRO_drugTargets.Phase3_N141I_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)

#Show the .svg file
showSVG(myNV_CADRO_drugTargets.Phase3, opacity=0.3,outFile = "/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_CERNO.CADRO_drugTargets.Phase3_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))


######################################################################
#49. Create Chord Diagram for subset of Phase 3 drug agents with defined
#CADRO mechanism classes and pathways across all mutations
######################################################################

#Filter by pathway and drug trial phase and then get unique drug agents
#to create matrix of 0's and 1's for input matrix

#Biological Adhesion
uniqueDrugPathway.biologicalAdhesion.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp biological adhesion") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.biologicalAdhesion.Phase3 <- data.frame(unique(uniqueDrugPathway.biologicalAdhesion.Phase3$Agent))

#Cell Fate Specification
uniqueDrugPathway.cellFate.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp cell fate commitment") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.cellFate.Phase3 <- data.frame(unique(uniqueDrugPathway.cellFate.Phase3$Agent))

#CNS Neuron Differentiation
uniqueDrugPathway.CNS_neuronDiff.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp central nervous system neuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.CNS_neuronDiff.Phase3 <- data.frame(unique(uniqueDrugPathway.CNS_neuronDiff.Phase3$Agent))

#Cerebral Cortex GABAergic Interneuron Differentiation
uniqueDrugPathway.cortex_GABAergic.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp cerebral cortex gabaergic interneuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.cortex_GABAergic.Phase3 <- data.frame(unique(uniqueDrugPathway.cortex_GABAergic.Phase3$Agent))

#G protein coupled receptor signaling
uniqueDrugPathway.Gprotein_signaling.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp g protein coupled receptor signaling pathway") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.Gprotein_signaling.Phase3 <- data.frame(unique(uniqueDrugPathway.Gprotein_signaling.Phase3$Agent))

#Nervous system process
uniqueDrugPathway.nervousSystem.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp nervous system process") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.nervousSystem.Phase3 <- data.frame(unique(uniqueDrugPathway.nervousSystem.Phase3$Agent))

#Neuron Fate Committment
uniqueDrugPathway.neuronFateCommit.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp neuron fate commitment") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.neuronFateCommit.Phase3 <- data.frame(unique(uniqueDrugPathway.neuronFateCommit.Phase3$Agent))

#Taxis
uniqueDrugPathway.taxis.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp taxis") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.taxis.Phase3 <- data.frame(unique(uniqueDrugPathway.taxis.Phase3$Agent))

#Neuron Fate Specification
uniqueDrugPathway.neuronFateSpec.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp neuron fate specification") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.neuronFateSpec.Phase3 <- data.frame(unique(uniqueDrugPathway.neuronFateSpec.Phase3$Agent))

#Regulation of neuron differentiation
uniqueDrugPathway.reg_neuronDiff.Phase3 <- chosenDrugPathways.CERNO.CADROclasses %>%
  dplyr::filter(Title == "Gobp regulation of neuron differentiation") %>%
  dplyr::filter(Phase_Trial == "Phase 3")
uniqueDrugAgent.reg_neuronDiff.Phase3 <- data.frame(unique(uniqueDrugPathway.reg_neuronDiff.Phase3$Agent))

#Read .csv file with drug input matrix 
chordDiagram_Values_drugTargets.Phase3 <- read.csv("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Common_Drug_Targets/CERNO.Pathways.CADRO.classes/230805_chordDiagValues_allMutations_drugTargets_Phase3.csv")

#Source: https://urldefense.com/v3/__https://stackoverflow.com/questions/73609132/how-to-create-a-chord-diagram-in-r__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqcD_PAlf$ 
#Set the colors for each drug
#Source: https://urldefense.com/v3/__https://r-charts.com/flow/chord-diagram/__;!!Mih3wA!GbFXsHnOACf0klwxladpINv8v89j8USAF-Bd1X2CnOpmeYryR13GTC76QkhPAH7lzPpxALW4XMcJGUagfRBNqe_Um70Q$ 
colors7.Phase3 <- c(Blarcamesine = "#cc6677",
                    Guanfacine = "#44aa99", 
                    Caffeine = "#88ccee",
                    Escitalopram = "#2C2E9C",
                    AVP.786 = "#b6dbff",
                    AXS.05 = "#FFC20A",
                    Octohydroaminoacridine.Succinate = "#FFB000",
                    Nabilone = "#FFC107",
                    Brexpiprazole = "#000000",
                    Donepezil = "#E69F00",
                    NE3107 = "#000000"
)

#Output to an .svg file and create the labelled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Phase3_Chord_labeled_v2.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_drugTargets.Phase3, grid.col = colors7.Phase3)
dev.off()

#Output to an .svg file and create the unlabeled chord diagram
svglite("/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/Drug_Target_Analysis/Chord_Diagrams/All3_PSEN1_2_APP_dATAC_DEGs_Drug_Targets_Phase3_Chord_unlabeled_v2.svg")
set.seed(123)
chordDiagram(chordDiagram_Values_drugTargets.Phase3, grid.col = colors7.Phase3, annotationTrack = "grid")
dev.off()


#############################   
#50. Save the .RData File
#############################

#Save the .RData File
save.image('/Users/phoebevaldes/Desktop/EOFAD_MS_Revisions/RData/230818_ATAC_tmodCERNO_drugTargets_Env.RData')
