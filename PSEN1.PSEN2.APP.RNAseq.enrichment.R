#RNA-Seq Enrichment Analysis pipeline 
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
library("fgsea")
library("tmod")
library("topconfects")
library("msigdbr")
library("viper")
library("dorothea")
# load limma-voom efit object
load("~/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")

############################################################
#Loading the Gene Set Databases
############################################################

#Load the geneset database - down loadable for GOBP, Hallmark, etc
#You already imported the Gene Ontology - Biological Process as Hs.GOBP
Hs.GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
Hs.Hallmark <- msigdbr(species = "Homo sapiens", category = "H")

#We need to select Symbol now instead of ENTREZ ID
Hs.GOBP.Symbol <- Hs.GOBP %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.GOBP.Symbol %>% 
  head() %>% 
  lapply(head)

#Select Hallmark Symbol
Hs.Hallmark.Symbol <- Hs.Hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Hallmark.Symbol %>% 
  head() %>% 
  lapply(head)

#Using your own gmt file; provide the gmt and filepath
Hs.ECC <-  qusage::read.gmt(file.path("~/gmtfiles", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt"))
Hs.ReMap <- qusage::read.gmt(file.path("~/gmtfiles", "ReMap_ChIP-seq_expanded.gmt"))
Hs.miR <- qusage::read.gmt(file.path("~/gmtfiles", "miRTarBase_2017.gmt"))

##################################################################
#T-value Limma Ranking for APP-V717I vs. NDC
##################################################################

# You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

#Convert your limma topTable results to a data table format
reslimma_V717I_dt <- data.table(reslimma_V717I)


#Remove rows with an empty or duplicated Symbol
reslimma_V717I_dt <- subset(reslimma_V717I_dt, SYMBOL != "" )
reslimma_V717I_dt <- subset(reslimma_V717I_dt, ! duplicated(SYMBOL))

#Tidy up and order by t metric
ranks_limma_V717I <- as_tibble(reslimma_V717I_dt[order(t), list(SYMBOL, t)])
ranks_limma_V717I

#Deframe the data table
ranks_t_V717I <- deframe(ranks_limma)

#see the top 20 genes by t
head(ranks_t_V717I, 20)

##################################################################
#T-value Limma Ranking for PSEN1-A79V vs. NDC
##################################################################

# You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Convert your limma topTable results to a data table format
reslimma_A79V_dt <- data.table(reslimma_A79V)


#Remove rows with an empty or duplicated Symbol
reslimma_A79V_dt <- subset(reslimma_A79V_dt, SYMBOL != "" )
reslimma_A79V_dt <- subset(reslimma_A79V_dt, ! duplicated(SYMBOL))

#Tidy up and order by t metric
ranks_limma_A79V <- as_tibble(reslimma_A79V_dt[order(t), list(SYMBOL, t)])
ranks_limma_A79V

#Deframe the data table
ranks_t_A79V <- deframe(ranks_limma_A79V)

#see the top 20 genes by t
head(ranks_t_A79V, 20)

##################################################################
#T-value Limma Ranking for PSEN2-N141I vs. NDC
##################################################################

# You could also use topTable at a particular FDR (those which dont have a logFC)
reslimma_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Convert your limma topTable results to a data table format
reslimma_N141I_dt <- data.table(reslimma_N141I)


#Remove rows with an empty or duplicated Symbol
reslimma_N141I_dt <- subset(reslimma_N141I_dt, SYMBOL != "" )
reslimma_N141I_dt <- subset(reslimma_N141I_dt, ! duplicated(SYMBOL))

#Tidy up and order by t metric
ranks_limma_N141I <- as_tibble(reslimma_N141I_dt[order(t), list(SYMBOL, t)])
ranks_limma_N141I

#Deframe the data table
ranks_t_N141I <- deframe(ranks_limma_N141I)

#see the top 20 genes by t
head(ranks_t_N141I, 20)

##################################################################################
#Running CERNO for efit object using t-value limma for GOBP/Hallmark 
##################################################################################

#Run enrichment using the CERNO test in the tmod package
#Coincident Extreme Ranks in Numerical Observations (CERNO)
#See Zyla et al. 2019 - https://academic.oup.com/bioinformatics/article/35/24/5146/5511403
#Ranked enrichment approach
#This allows for simultaneous enrichment of all contrasts in your limma contrast matrix
#Does not differentiate between up- or down-regulated pathways, but can add DEG status of all genes within enriched pathway


#Import the latest MSigDB XML File -> download 7.4 version
msig <- tmodImportMSigDB(file.path("~/msigdb", "msigdb_v7.4.xml"))

#Select the GO BP gene sets (you can select whichever geneset db you would like)
unique(as.character(msig$MODULES$Subcategory))
unique(as.character(msig$MODULES$Category))
Hallmark.sel <- msig$MODULES$Category == "H"
GOBP.sel <- msig$MODULES$Subcategory == "C5_GO:BP"
Reactome.sel <- msig$MODULES$Subcategory == "CP:REACTOME"
GOBP.exp.sel <- msig$MODULES$Subcategory == "GO:BP"

#Run the CERNO enrichment test on a limma object which you've applied lmfit and eBayes
#Supply the gene annotation - here, which is contained in the voom object v - v$genes$SYMBOL
#Sort by Minimum Significant Distance (msd)
#MSD is the signed distance of the confidence interval (CI) of the logarithm of fold change (logFC) estimate from no change (zero)

#filter pathways on certain adjusted p-value -> FDR
#use dplyr and subset

CERNO.GOBP <- tmodLimmaTest(efit, 
                            v$genes$SYMBOL, 
                            sort.by = "msd", 
                            tmodFunc =tmodCERNOtest,
                            mset=msig[GOBP.sel])

CERNO.Hallmark <- tmodLimmaTest(efit, 
                                v$genes$SYMBOL, 
                                sort.by = "msd", 
                                tmodFunc =tmodCERNOtest,
                                mset=msig[Hallmark.sel])

#See the top results for all contrasts within your study
head(CERNO.GOBP)
head(CERNO.Hallmark)

write.csv(CERNO.GOBP$V717IvsNDC, file = "~/CERNO_GOBP_Results_V717IvsNDC_eBayes.csv") 
write.csv(CERNO.GOBP$A79VvsNDC, file = "~/CERNO_GOBP_Results_A79VvsNDC_eBayes.csv") 
write.csv(CERNO.GOBP$N141IvsNDC, file = "~/CERNO_GOBP_Results_N141IvsNDC_eBayes.csv") 
write.csv(CERNO.Hallmark$V717IvsNDC, file = "~/CERNO_Hallmark_Results_V717IvsNDC_eBayes.csv") 
write.csv(CERNO.Hallmark$A79VvsNDC, file = "~/CERNO_Hallmark_Results_A79VvsNDC_eBayes.csv") 
write.csv(CERNO.Hallmark$N141IvsNDC, file = "~/CERNO_Hallmark_Results_N141IvsNDC_eBayes.csv") 

##################################################################################
#Running CERNO for efit object using t-value limma for ENCODE/ChEA 
##################################################################################

#------------------------------------------------------------------------#
# Make ENCODE/ChEA Consensus mset for use with tmod 
ECC <- read_table("~/gmtfiles/ECC_TFchipenrich.txt")
ECC <- ECC[,c(2,1)]
colnames(ECC) <- c("gene_id","go_id")
m2g_ECC <- split(ECC$gene_id, ECC$go_id)
#gt <- toTable(GOTERM)
m_ECC <- data.frame(ID=names(m2g_ECC))
m_ECC$Title <- m_ECC$ID
goset_ECC <- makeTmod(modules=m_ECC, modules2genes=m2g_ECC)
msetECC <- goset_ECC
save(msetECC,file="~/gmtfiles/msetECC.Rdata")

# Run CERNO
CERNO.ECC <- tmodLimmaTest(efit,v_new$genes$SYMBOL, sort.by = "msd", tmodFunc = tmodCERNOtest, mset= msetECC)
head(CERNO.ECC)

#Save the results using the eBayes object
write.csv(CERNO.ECC$V717IvsNDC, file = "~/CERNO_ECC_Results_V717IvsNDC_eBayes.csv") 
write.csv(CERNO.ECC$A79VvsNDC, file = "~/CERNO_ECC_Results_A79VvsNDC_eBayes.csv") 
write.csv(CERNO.ECC$N141IvsNDC, file = "~/CERNO_ECC_Results_N141IvsNDC_eBayes.csv") 

##################################################################################
#Running CERNO for efit object using t-value limma for ReMap 
##################################################################################

#------------------------------------------------------------------------#
# Make ReMap ChIP-Seq mset for use with tmod 

#Note: make sure to compare ReMap.MODULES2GENES and ReMap.MODULES first

#Load ReMap.gmt files
ReMap.gmt <- read_table("~/gmtfiles/ReMap_chipenrich.txt")
ReMap <- ReMap[,c(2,1)]
colnames(ReMap) <- c("gene_id","go_id")
m2g_ReMap <- split(ReMap$gene_id, ReMap$go_id)
#gt <- toTable(GOTERM)
m_ReMap <- data.frame(ID=names(m2g_ReMap))
m_ReMap$Title <- m_ReMap$ID
goset_ReMap <- makeTmod(modules=m_ReMap, modules2genes=m2g_ReMap)
msetReMap <- goset_ReMap
save(msetReMap,file="~/gmtfiles/msetReMap.Rdata")

# Run CERNO
CERNO.ReMap <- tmodLimmaTest(efit,v_new$genes$SYMBOL, sort.by = "msd", tmodFunc = tmodCERNOtest, mset=msetReMap)
head(CERNO.ReMap)

#Save the results using eBayes object
write.csv(CERNO.ReMap$V717IvsNDC, file = "~/CERNO_ReMap_Results_V717IvsNDC_eBayes.csv") 
write.csv(CERNO.ReMap$A79VvsNDC, file = "~/CERNO_ReMap_Results_A79VvsNDC_eBayes.csv") 
write.csv(CERNO.ReMap$N141IvsNDC, file = "~/CERNO_ReMap_Results_N141IvsNDC_eBayes.csv") 

################################################################
#fGSEA with t-value Limma Ranking for GOBP for APP-V717I vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_V717I as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes_V717I_t <- fgsea::fgseaMultilevel(pathways=Hs.GOBP.Symbol, stats=ranks_t_V717I, absEps = 0)

# Make the Results tidy
fgseaResTidy_V717I_t <- fgseaRes_V717I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy_V717I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy_V717I_t <- fgseaResTidy_V717I_t[order(fgseaResTidy_V717I_t$padj),]
fgseaResTidy_V717I_t_dt <- data.table(fgseaResTidy_V717I_t)
fgseaResTidy_V717I_t_df <- apply(fgseaResTidy_V717I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy_V717I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_GOBP_V717I_NDC_t.csv", row.names = TRUE)

################################################################
#fGSEA with t-value Limma Ranking for GOBP for PSEN1-A79V vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_A79V as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes_A79V_t <- fgsea::fgseaMultilevel(pathways=Hs.GOBP.Symbol, stats=ranks_t_A79V, absEps = 0)

# Make the Results tidy
fgseaResTidy_A79V_t <- fgseaRes_A79V_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy_A79V_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy_A79V_t <- fgseaResTidy_A79V_t[order(fgseaResTidy_A79V_t$padj),]
fgseaResTidy_A79V_t_dt <- data.table(fgseaResTidy_A79V_t)
fgseaResTidy_A79V_t_df <- apply(fgseaResTidy_A79V_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy_A79V_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_GOBP_A79V_NDC_t.csv", row.names = TRUE)

################################################################
#fGSEA with t-value Limma Ranking for GOBP for PSEN2-N141I vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_N141I as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes_N141I_t <- fgsea::fgseaMultilevel(pathways=Hs.GOBP.Symbol, stats=ranks_t_N141I, absEps = 0)

# Make the Results tidy
fgseaResTidy_N141I_t <- fgseaRes_N141I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy_N141I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy_N141I_t <- fgseaResTidy_N141I_t[order(fgseaResTidy_N141I_t$padj),]
fgseaResTidy_N141I_t_dt <- data.table(fgseaResTidy_N141I_t)
fgseaResTidy_N141I_t_df <- apply(fgseaResTidy_N141I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy_N141I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_GOBP_N141I_NDC_t.csv", row.names = TRUE)

#####################################################################
#fGSEA with t-value Limma Ranking for Hallmark for APP-V717I vs. NDC
#####################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_V717I as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes.Hallmark_V717I_t <- fgsea::fgseaMultilevel(pathways=Hs.Hallmark.Symbol, stats=ranks_t_V717I, absEps = 0)

# Make the Results tidy
fgseaResTidy.Hallmark_V717I_t <- fgseaRes.Hallmark_V717I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.Hallmark_V717I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.Hallmark_V717I_t <- fgseaResTidy.Hallmark_V717I_t[order(fgseaResTidy.Hallmark_V717I_t$padj),]
fgseaResTidy.Hallmark_V717I_t_dt <- data.table(fgseaResTidy.Hallmark_V717I_t)
fgseaResTidy.Hallmark_V717I_t_df <- apply(fgseaResTidy.Hallmark_V717I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.Hallmark_V717I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_Hallmark_V717I_NDC_t.csv", row.names = TRUE)

#####################################################################
#fGSEA with t-value Limma Ranking for Hallmark for PSEN1-A79V vs. NDC
#####################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_A79V as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes.Hallmark_A79V_t <- fgsea::fgseaMultilevel(pathways=Hs.Hallmark.Symbol, stats=ranks_t_A79V, absEps = 0)

# Make the Results tidy
fgseaResTidy.Hallmark_A79V_t <- fgseaRes.Hallmark_A79V_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.Hallmark_A79V_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.Hallmark_A79V_t <- fgseaResTidy.Hallmark_A79V_t[order(fgseaResTidy.Hallmark_A79V_t$padj),]
fgseaResTidy.Hallmark_A79V_t_dt <- data.table(fgseaResTidy.Hallmark_A79V_t)
fgseaResTidy.Hallmark_A79V_t_df <- apply(fgseaResTidy.Hallmark_A79V_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.Hallmark_A79V_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_Hallmark_A79V_NDC_t.csv", row.names = TRUE)

#####################################################################
#fGSEA with t-value Limma Ranking for Hallmark for PSEN2-N141I vs. NDC
#####################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_N141I as a ranking for fGSEA 

# Run multilevel fgsea (more accurate p-value determination - takes longer but better)
fgseaRes.Hallmark_N141I_t <- fgsea::fgseaMultilevel(pathways=Hs.Hallmark.Symbol, stats=ranks_t_N141I, absEps = 0)

# Make the Results tidy
fgseaResTidy.Hallmark_N141I_t <- fgseaRes.Hallmark_N141I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.Hallmark_N141I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.Hallmark_N141I_t <- fgseaResTidy.Hallmark_N141I_t[order(fgseaResTidy.Hallmark_N141I_t$padj),]
fgseaResTidy.Hallmark_N141I_t_dt <- data.table(fgseaResTidy.Hallmark_N141I_t)
fgseaResTidy.Hallmark_N141I_t_df <- apply(fgseaResTidy.Hallmark_N141I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.Hallmark_N141I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_Hallmark_N141I_NDC_t.csv", row.names = TRUE)

################################################################
#fGSEA with t-value Limma Ranking for ECC for APP-V717I vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_V717I as a ranking for fGSEA 

#Run multilevel fgsea (more accurate p-value determination - takes longer but better)
#Run ECC enrichment for using fgsea (make sure to set a nproc less than the total processors you have available)
fgseaRes.ECC_V717I_t <- fgsea::fgseaMultilevel(pathways=Hs.ECC, stats=ranks_t_V717I, absEps = 0, nproc = 12)

# Make the Results tidy
fgseaResTidy.ECC_V717I_t <- fgseaRes.ECC_V717I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.ECC_V717I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.ECC_V717I_t <- fgseaResTidy.ECC_V717I_t[order(fgseaResTidy.ECC_V717I_t$padj),]
fgseaResTidy.ECC_V717I_t_dt <- data.table(fgseaResTidy.ECC_V717I_t)
fgseaResTidy.ECC_V717I_t_df <- apply(fgseaResTidy.ECC_V717I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.ECC_V717I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_ECC_V717I_NDC_t_v3.csv", row.names = TRUE)

################################################################
#fGSEA with t-value Limma Ranking for ECC for PSEN1-A79V vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_A79V as a ranking for fGSEA 

#Run multilevel fgsea (more accurate p-value determination - takes longer but better)
#Run ECC enrichment for using fgsea (make sure to set a nproc less than the total processors you have available)
fgseaRes.ECC_A79V_t <- fgsea::fgseaMultilevel(pathways=Hs.ECC, stats=ranks_t_A79V, absEps = 0, nproc = 12)

# Make the Results tidy
fgseaResTidy.ECC_A79V_t <- fgseaRes.ECC_A79V_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.ECC_A79V_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.ECC_A79V_t <- fgseaResTidy.ECC_A79V_t[order(fgseaResTidy.ECC_A79V_t$padj),]
fgseaResTidy.ECC_A79V_t_dt <- data.table(fgseaResTidy.ECC_A79V_t)
fgseaResTidy.ECC_A79V_t_df <- apply(fgseaResTidy.ECC_A79V_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.ECC_A79V_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_ECC_A79V_NDC_t_v3.csv", row.names = TRUE)

################################################################
#fGSEA with t-value Limma Ranking for ECC for PSEN2-N141I vs. NDC
################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_A79V as a ranking for fGSEA 

#Run multilevel fgsea (more accurate p-value determination - takes longer but better)
#Run ECC enrichment for using fgsea (make sure to set a nproc less than the total processors you have available)
fgseaRes.ECC_N141I_t <- fgsea::fgseaMultilevel(pathways=Hs.ECC, stats=ranks_t_N141I, absEps = 0, nproc = 12)

# Make the Results tidy
fgseaResTidy.ECC_N141I_t <- fgseaRes.ECC_N141I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.ECC_N141I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.ECC_N141I_t <- fgseaResTidy.ECC_N141I_t[order(fgseaResTidy.ECC_N141I_t$padj),]
fgseaResTidy.ECC_N141I_t_dt <- data.table(fgseaResTidy.ECC_N141I_t)
fgseaResTidy.ECC_N141I_t_df <- apply(fgseaResTidy.ECC_N141I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.ECC_N141I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_ECC_N141I_NDC_t_v3.csv", row.names = TRUE)

#################################################################
#fGSEA with t-value Limma Ranking for ReMap for APP-V717I vs. NDC
#################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_V717I as a ranking for fGSEA 

#Run multilevel fgsea (more accurate p-value determination - takes longer but better)
#Run ECC enrichment for using fgsea (make sure to set a nproc less than the total processors you have available)
fgseaRes.ReMap_V717I_t <- fgsea::fgseaMultilevel(pathways=Hs.ReMap, stats=ranks_t_V717I, absEps = 0, nproc = 12)

# Make the Results tidy
fgseaResTidy.ReMap_V717I_t <- fgseaRes.ReMap_V717I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.ReMap_V717I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.ReMap_V717I_t <- fgseaResTidy.ReMap_V717I_t[order(fgseaResTidy.ReMap_V717I_t$padj),]
fgseaResTidy.ReMap_V717I_t_dt <- data.table(fgseaResTidy.ReMap_V717I_t)
fgseaResTidy.ReMap_V717I_t_df <- apply(fgseaResTidy.ReMap_V717I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.ReMap_V717I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_ReMap_V717I_NDC_t.csv", row.names = TRUE)

##################################################################
#fGSEA with t-value Limma Ranking for ReMap for PSEN1-A79V vs. NDC
##################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_A79V as a ranking for fGSEA

#Run multilevel fgsea (more accurate p-value determination - takes longer but better)
#Run ECC enrichment for using fgsea (make sure to set a nproc less than the total processors you have available)
fgseaRes.ReMap_A79V_t <- fgsea::fgseaMultilevel(pathways=Hs.ReMap, stats=ranks_t_A79V, absEps = 0, nproc = 12)

# Make the Results tidy
fgseaResTidy.ReMap_A79V_t <- fgseaRes.ReMap_A79V_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.ReMap_A79V_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.ReMap_A79V_t <- fgseaResTidy.ReMap_A79V_t[order(fgseaResTidy.ReMap_A79V_t$padj),]
fgseaResTidy.ReMap_A79V_t_dt <- data.table(fgseaResTidy.ReMap_A79V_t)
fgseaResTidy.ReMap_A79V_t_df <- apply(fgseaResTidy.ReMap_A79V_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.ReMap_A79V_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_ReMap_A79V_NDC_t.csv", row.names = TRUE)

###################################################################
#fGSEA with t-value Limma Ranking for ReMap for PSEN2-N141I vs. NDC
###################################################################

#--------------------------------------------------------------------------------------------------------
#Use ranks_t_N141I as a ranking for fGSEA 

#Run multilevel fgsea (more accurate p-value determination - takes longer but better)
#Run ECC enrichment for using fgsea (make sure to set a nproc less than the total processors you have available)
fgseaRes.ReMap_N141I_t <- fgsea::fgseaMultilevel(pathways=Hs.ReMap, stats=ranks_t_N141I, absEps = 0, nproc = 12)

# Make the Results tidy
fgseaResTidy.ReMap_N141I_t <- fgseaRes.ReMap_N141I_t %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table
fgseaResTidy.ReMap_N141I_t %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Arrange table from smallest to largest adjusted p-value
fgseaResTidy.ReMap_N141I_t <- fgseaResTidy.ReMap_N141I_t[order(fgseaResTidy.ReMap_N141I_t$padj),]
fgseaResTidy.ReMap_N141I_t_dt <- data.table(fgseaResTidy.ReMap_N141I_t)
fgseaResTidy.ReMap_N141I_t_df <- apply(fgseaResTidy.ReMap_N141I_t_dt,2,as.character)

dev.off()

#Save results to .csv file
write.csv(fgseaResTidy.ReMap_N141I_t_df, file = "~/fGSEA_t-value_limma/Results_fGSEA_ReMap_N141I_NDC_t.csv", row.names = TRUE)

####################################################################
#Mitch rank-MANOVA enrichment analysis
####################################################################

Hs.GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
GOBP.genesets <- Hs.GOBP %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Hallmark <- msigdbr(species = "Homo sapiens", category = "H")
Hallmark.genesets <- Hs.Hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
#
# Run topconfects
confects.RNA.A79V.filter  <- limma_confects(efit, coef="A79VvNDC", fdr=0.05)
confects.RNA.V717I.filter <- limma_confects(efit, coef="V717IvNDC", fdr=0.05)
confects.RNA.N141I.filter <- limma_confects(efit, coef="N141IvNDC", fdr=0.05)
# remove NAs
confects.RNA.V717I.symbol <- na.omit(confects.RNA.V717I.filter[,c(2:3)])
confects.RNA.A79V.symbol <- na.omit(confects.RNA.A79V.filter[,c(2:3)])
confects.RNA.N141I.symbol <- na.omit(confects.RNA.N141I.filter[,c(2:3)])
#remove duplicates
confects.RNA.V717I.symbol <- confects.RNA.V717I.symbol[!duplicated(confects.RNA.V717I.symbol$SYMBOL),]
confects.RNA.A79V.symbol <- confects.RNA.A79V.symbol[!duplicated(confects.RNA.A79V.symbol$SYMBOL),]
confects.RNA.N141I.symbol <- confects.RNA.N141I.symbol[!duplicated(confects.RNA.N141I.symbol$SYMBOL),]
# Reneame rows and columns
rownames(confects.RNA.V717I.symbol) <- confects.RNA.V717I.symbol$SYMBOL
rownames(confects.RNA.A79V.symbol) <- confects.RNA.A79V.symbol$SYMBOL
rownames(confects.RNA.N141I.symbol) <- confects.RNA.N141I.symbol$SYMBOL
colnames(confects.RNA.V717I.symbol)[2] <- "confect"
colnames(confects.RNA.A79V.symbol)[2] <- "confect"
colnames(confects.RNA.N141I.symbol)[2] <- "confect"
# Make mitch enrichment dataframe
x.mitch.RNA <-list("V717I"=as.data.frame(confects.RNA.V717I.symbol[,c(2)]),
                   "A79V"=as.data.frame(confects.RNA.A79V.symbol[,c(2)]),
                   "N141I"=as.data.frame(confects.RNA.N141I.symbol[,c(2)]))
rownames(x.mitch.RNA$V717I) <- confects.RNA.V717I.symbol$SYMBOL
rownames(x.mitch.RNA$A79V) <- confects.RNA.A79V.symbol$SYMBOL
rownames(x.mitch.RNA$N141I) <- confects.RNA.N141I.symbol$SYMBOL
colnames(x.mitch.RNA$V717I) <- "confect"
colnames(x.mitch.RNA$A79V) <- "confect"
colnames(x.mitch.RNA$N141I) <- "confect"
# Import topconfects prescored rankings
y.mitch.RNA <-mitch_import(x.mitch.RNA,DEtype="prescored")
# Run Hallmark enrichment
all3.RNA <-mitch_calc(y.mitch.RNA,
                      Hallmark.genesets,
                      priority="significance",
                      cores=4)
write.csv(all3.RNA$enrichment_result,file="PSEN1_PSEN2_APP_topconfects_mitch_Hallmark.csv")
# Run GOBP enrichment
all3.RNA.GOBP <-mitch_calc(y.mitch.RNA,
                           GOBP.genesets,
                           priority="significance",
                           cores=4)
write.csv(all3.RNA.GOBP$enrichment_result,file="PSEN1_PSEN2_APP_topconfects_mitch_GOBP.csv")

####################################################################
#Dorthea analysis for APP-V717I vs.NDC
####################################################################

#---------------------------------------------------------------------------------------------------
## DoRothEA TF Activity Analysis ## 
# What if you want to determine the activity of a given TF? Use DoRothEA
# TF regulon genesets have been pre-generated by Garcia-Alonso et al 2019 Genome Research, 29(8), 1363-1375
# https://doi.org/10.1101/gr.240663.118
# See https://github.com/saezlab/DoRothEA

# Load TF regulon genesets in VIPER format
# There are 6 sets A-E (plus the top 10) in descending stringency (E being the least stringent, with computational predictions)
# Download from github, or from $RNASEQ
# C set is suitable - it requires either 1) Curated data w/ TFBS or 2) ChIP-Seq data w/ TFBS
load('~/DoRothEA_TFs/C_viperRegulon.rdata')

# Clean TF names & explore object
names(viper_regulon) <- sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)

# Explore the regulons object
names(viper_regulon)[1:10]
viper_regulon[[1]]

#Perform APP-V717I vs. NDC Analysis
# This assumes you have done limma or limma-voom analysis for microarray or RNA-Seq. You can customize if you use DESeq2, EdgeR, etc
# You can choose to order by z-scored (sign(logFC) x p.val)

# Create topTable results 
res_limma_dr_V717I <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf)

#Use limma topTable for ranking
DEsignature_V717I <- res_limma_dr_V717I

#Match gene symbols to Ensembl ID's
DEsignature_V717I = merge(DEsignature_V717I, Hs_ann, by="SYMBOL")

#Convert your limma topTable results to a data table format
DEsignature_V717I <- data.table(DEsignature_V717I)

# Exclude probes with unknown or duplicated gene symbol
DEsignature_V717I <- subset(DEsignature_V717I, SYMBOL != "" )
DEsignature_V717I <- subset(DEsignature_V717I, ! duplicated(SYMBOL))

# Estimate z-score values for the GES. Check VIPER manual for details
myStatistics_V717I <- matrix(DEsignature_V717I$logFC, dimnames = list(DEsignature_V717I$SYMBOL, 'logFC') )

myQvalue_V717I <- matrix(DEsignature_V717I$adj.P.Val, dimnames = list(DEsignature_V717I$SYMBOL, 'adj.P.Val') )
mySignature_V717I <- (qnorm(myQvalue_V717I/2, lower.tail = FALSE) * sign(myStatistics_V717I))[, 1]

# Reorder and rename
limmaSignature_V717I <- mySignature_V717I[order(mySignature_V717I, decreasing = T)]

# Estimate TF activities
# Use limmaSignature
mrs_V717I <- msviper(ges = limmaSignature_V717I, regulon = viper_regulon, minsize = 4, ges.filter = F, verbose = TRUE)

# Find the top contributing genes to each regulon
mrs_V717I <- ledge(mrs_V717I)

# See the results
summary(mrs_V717I)

# Save the results as a dataframe
TF_activities_V717I <- data.frame(Regulon = names(mrs_V717I$es$nes),
                                  Size = mrs_V717I$es$size[ names(mrs_V717I$es$nes) ], 
                                  NES = mrs_V717I$es$nes, 
                                  p.value = mrs_V717I$es$p.value, 
                                  FDR = p.adjust(mrs_V717I$es$p.value, method = 'fdr'))
# Tidy up and order by FDR
TF_activitiesTidy_V717I <- TF_activities_V717I %>%
  as_tibble() %>%
  arrange(FDR)

# Show in a nice table
DT::datatable(TF_activitiesTidy_V717I)

#Arrange table from smallest to largest adjusted p-value
TF_activitiesTidy_V717I <- TF_activitiesTidy_V717I[order(TF_activitiesTidy_V717I$FDR),]

# Save results
write.csv(TF_activitiesTidy_V717I, file = '~/DoRothEAvC_TF_activities_eBayes_V717I_NDC_qValue.csv')

####################################################################
#Performing Dorthea analysis for PSEN1-A79V vs.NDC
####################################################################

#---------------------------------------------------------------------------------------------------
#Perform PSEN1-A79V vs. NDC Analysis
# This assumes you have done limma or limma-voom analysis for microarray or RNA-Seq. You can customize if you use DESeq2, EdgeR, etc
# You can choose to order by z-scored (sign(logFC) x p.val)

# Create topTable results 
res_limma_dr_A79V <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf)

#Use limma topTable for ranking
DEsignature_A79V <- res_limma_dr_A79V

#Match gene symbols to Ensembl ID's
DEsignature_A79V = merge(DEsignature_A79V, Hs_ann, by="SYMBOL")

#Convert your limma topTable results to a data table format
DEsignature_A79V <- data.table(DEsignature_A79V)

# Exclude probes with unknown or duplicated gene symbol
DEsignature_A79V <- subset(DEsignature_A79V, SYMBOL != "" )
DEsignature_A79V <- subset(DEsignature_A79V, ! duplicated(SYMBOL))

# Estimate z-score values for the GES. Check VIPER manual for details
myStatistics_A79V <- matrix(DEsignature_A79V$logFC, dimnames = list(DEsignature_A79V$SYMBOL, 'logFC') )
myQvalue_A79V <- matrix(DEsignature_A79V$adj.P.Val, dimnames = list(DEsignature_A79V$SYMBOL, 'adj.P.Val') )
mySignature_A79V <- (qnorm(myQvalue_A79V/2, lower.tail = FALSE) * sign(myStatistics_A79V))[, 1]

# Reorder and rename
limmaSignature_A79V <- mySignature_A79V[order(mySignature_A79V, decreasing = T)]

# Estimate TF activities
# Use limmaSignature
mrs_A79V <- msviper(ges = limmaSignature_A79V, regulon = viper_regulon, minsize = 4, ges.filter = F, verbose = TRUE)

# Find the top contributing genes to each regulon
mrs_A79V <- ledge(mrs_A79V)

# See the results
summary(mrs_A79V)

# Save the results as a dataframe
TF_activities_A79V <- data.frame(Regulon = names(mrs_A79V$es$nes),
                                 Size = mrs_A79V$es$size[ names(mrs_A79V$es$nes) ], 
                                 NES = mrs_A79V$es$nes, 
                                 p.value = mrs_A79V$es$p.value, 
                                 FDR = p.adjust(mrs_A79V$es$p.value, method = 'fdr'))
# Tidy up and order by FDR
TF_activitiesTidy_A79V <- TF_activities_A79V %>%
  as_tibble() %>%
  arrange(FDR)

# Show in a nice table
DT::datatable(TF_activitiesTidy_A79V)

#Arrange table from smallest to largest adjusted p-value
TF_activitiesTidy_A79V <- TF_activitiesTidy_A79V[order(TF_activitiesTidy_A79V$FDR),]

# Save results
write.csv(TF_activitiesTidy_A79V, file = '~/DoRothEAvC_TF_activities_eBayes_A79V_NDC_Qvalue.csv')

####################################################################
#Performing Dorothea analysis for PSEN2-N141I vs.NDC
####################################################################

#---------------------------------------------------------------------------------------------------
#Perform PSEN2-N141I vs. NDC Analysis
# This assumes you have done limma or limma-voom analysis for microarray or RNA-Seq. You can customize if you use DESeq2, EdgeR, etc
# You can choose to order by z-scored (sign(logFC) x p.val)
# Create topTable results 

res_limma_dr_N141I <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf)

#Use limma topTable for ranking
DEsignature_N141I <- res_limma_dr_N141I

#Match gene symbols to Ensembl ID's
DEsignature_N141I = merge(DEsignature_N141I, Hs_ann, by="SYMBOL")

#Convert your limma topTable results to a data table format
DEsignature_N141I <- data.table(DEsignature_N141I)

# Exclude probes with unknown or duplicated gene symbol
DEsignature_N141I <- subset(DEsignature_N141I, SYMBOL != "" )
DEsignature_N141I <- subset(DEsignature_N141I, ! duplicated(SYMBOL))

# Estimate z-score values for the GES. Check VIPER manual for details
myStatistics_N141I <- matrix(DEsignature_N141I$logFC, dimnames = list(DEsignature_N141I$SYMBOL, 'logFC') )
myQvalue_N141I <- matrix(DEsignature_N141I$adj.P.Val, dimnames = list(DEsignature_N141I$SYMBOL, 'adj.P.Val') )
mySignature_N141I <- (qnorm(myQvalue_N141I/2, lower.tail = FALSE) * sign(myStatistics_N141I))[, 1]

# Reorder and rename
limmaSignature_N141I <- mySignature_N141I[order(mySignature_N141I, decreasing = T)]

# Estimate TF activities
# Use limmaSignature
mrs_N141I <- msviper(ges = limmaSignature_N141I, regulon = viper_regulon, minsize = 4, ges.filter = F, verbose = TRUE)

# Find the top contributing genes to each regulon
mrs_N141I <- ledge(mrs_N141I)

# See the results
summary(mrs_N141I)

# Save the results as a dataframe
TF_activities_N141I <- data.frame(Regulon = names(mrs_N141I$es$nes),
                                  Size = mrs_N141I$es$size[ names(mrs_N141I$es$nes) ], 
                                  NES = mrs_N141I$es$nes, 
                                  p.value = mrs_N141I$es$p.value, 
                                  FDR = p.adjust(mrs_N141I$es$p.value, method = 'fdr'))
# Tidy up and order by FDR
TF_activitiesTidy_N141I <- TF_activities_N141I %>%
  as_tibble() %>%
  arrange(FDR)

# Show in a nice table
DT::datatable(TF_activitiesTidy_N141I)

#Arrange table from smallest to largest adjusted p-value
TF_activitiesTidy_N141I <- TF_activitiesTidy_N141I[order(TF_activitiesTidy_N141I$FDR),]

# Save results
write.csv(TF_activitiesTidy_N141I, file = '~/DoRothEAvC_TF_activities_eBayes_N141I_NDC_Qvalue.csv')