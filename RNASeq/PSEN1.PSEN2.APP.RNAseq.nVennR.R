#RNA-Seq nVenn figures
#Valdes et. al 2023 Molecular Psychiatry Submission
#Install the nVennR package
install.packages("nVennR")
#Load packages
library("limma")
library("VennDiagram")
library("nVennR")

#Load the efit object
load("~/PSEN1.PSEN2.APP.limma-voom.efit.Rdata")

#Do topTable() for the APP-V717I vs. NDC 
V717I_DEG <- topTable(efit, coef=1, adjust.method="BH", sort.by = "t", n = Inf) 
V717I_DEG <- V717I_DEG %>%
  dplyr::filter(adj.P.Val < 0.05)
V717I <- list(V717I = as.character(unique(V717I_DEG$SYMBOL)))

#Do topTable() for the PSEN1-A79V vs. NDC 
A79V_DEG <- topTable(efit, coef=2, adjust.method="BH", sort.by = "t", n = Inf) 
A79V_DEG <- A79V_DEG %>%
  dplyr::filter(adj.P.Val < 0.05)
A79V <- list(A79V = as.character(unique(A79V_DEG$SYMBOL)))

#Do topTable() for the PSEN2-N141I vs. NDC  
N141I_DEG <- topTable(efit, coef=3, adjust.method="BH", sort.by = "t", n = Inf) 
N141I_DEG <- N141I_DEG %>%
  dplyr::filter(adj.P.Val < 0.05)
N141I <- list(N141I = as.character(unique(N141I_DEG$SYMBOL)))

#Plot the nVenn Diagram
myNV3 <- plotVenn(list(V717I, A79V, N141I), sNames=c("V717I", "A79V", "N141I"),showPlot = T,nCycles = 1000000)
showSVG(myNV3, opacity=0.3,outFile = "~/PSEN1.PSEN2.APP.DEG_nVenn.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))