#RNA-seq RRHO2 analysis
#Valdes et. al 2023 Molecular Psychiatry Submission
#-------------------------------------------------------#
BiocManager::install("GeneOverlap")
library("devtools")
install_github("RRHO2/RRHO2")
library("GeneOverlap")
library("RRHO2")
# load efit object
load("~/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")
# Get DEGs for each comparison
V717I_DEGs <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t",p.value=0.05, n = Inf)
A79V_DEGs <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", p.value=0.05, n = Inf)
N141I_DEGs <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", p.value=0.05, n = Inf)
#
# test gene overlaps
go.obj.V717I.A79V <- newGeneOverlap(V717I_DEGs$ENSEMBL,
                                    A79V_DEGs$ENSEMBL,
                                    genome.size = length(res_limma_V717I$ENSEMBL))
go.obj.V717I.A79V <- testGeneOverlap(go.obj.V717I.A79V)
print(go.obj.V717I.A79V)
#
go.obj.N141I.A79V <- newGeneOverlap(N141I_DEGs$ENSEMBL,
                                    A79V_DEGs$ENSEMBL,
                                    genome.size = length(res_limma_N141I$ENSEMBL))
go.obj.N141I.A79V <- testGeneOverlap(go.obj.N141I.A79V)
print(go.obj.N141I.A79V)
#
go.obj.V717I.N141I <- newGeneOverlap(V717I_DEGs$ENSEMBL,
                                     N141I_DEGs$ENSEMBL,
                                     genome.size = length(res_limma_V717I$ENSEMBL))
go.obj.V717I.N141I <- testGeneOverlap(go.obj.V717I.N141I)
print(go.obj.V717I.N141I)
#
# RRHO2 overlap
# define new statistic
res_limma_V717I$stat <- ((-log10(res_limma_V717I$P.Value)) * sign(res_limma_V717I$logFC))
res_limma_A79V$stat <- ((-log10(res_limma_A79V$P.Value)) * sign(res_limma_A79V$logFC))
res_limma_N141I$stat <- ((-log10(res_limma_N141I$P.Value)) * sign(res_limma_N141I$logFC))
# Create list for each 
A79V_list <- res_limma_A79V[,c(1,10)]
V717I_list <- res_limma_V717I[,c(1,10)]
N141I_list <- res_limma_N141I[,c(1,10)]
# define viridis Gradient
viridisGradient <- viridis(101)
# Pairwise RRHO comparison of conditions
# PESN1 vs APP
RRHO_obj_A79V_V717I <- RRHO2::RRHO2_initialize(A79V_list,V717I_list, labels= c("A79V","V717I"), log10.ind = TRUE)
RRHO2_heatmap(RRHO_obj_A79V_V717I, colorGradient = viridisGradient)
RRHO_obj_A79V_V717I$hypermat
# PSEN2 vs APP
RRHO_obj_N141I_V717I <- RRHO2::RRHO2_initialize(N141I_list,V717I_list, labels= c("N141I","V717I"), log10.ind = TRUE)
dev.off()
RRHO2_heatmap(RRHO_obj_N141I_V717I, colorGradient = viridisGradient)
RRHO_obj_N141I_V717I$hypermat
# PSEN1 vs PSEN2
RRHO_obj_A79V_N141I <- RRHO2::RRHO2_initialize(A79V_list,N141I_list, labels= c("A79V","N141I"), log10.ind = TRUE)
dev.off()
RRHO2_heatmap(RRHO_obj_A79V_N141I, colorGradient = viridisGradient)
RRHO_obj_A79V_N141I$hypermat