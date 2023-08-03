#ATAC-seq/GWAS integration Analysis pipeline 
#Valdez et. al 2023 Molecular Psychiatry Submission
#Load packages
library("plyranges")
library("granges")
library("chipenrich")
library("ChIPseeker")

#load the data
load("~/ATAC-seq/dATACpeaks.V717I.annotated.Rdata")
load("~/ATAC-seq/dATACpeaks.A79V.annotated.Rdata")
load("~/ATAC-seq/dATACpeaks.N141I.annotated.Rdata")


####################################################################################
# Running the ChIPenrich program using R (V717I vs. NDC)
####################################################################################

#16,652 peaks found
dATACpeaks.V717I <- anno.V717I.filtered

#Load the differential ATAC peaks for APP-V717I

#16,652 peaks
length(dATACpeaks.V717I)

dATACpeaks.V717I.promoter <- dATACpeaks.V717I[dATACpeaks.V717I$annotation == "Promoter"]
#5,065 peaks
length(dATACpeaks.V717I.promoter)

dATACpeaks.V717I.nonpromoter <- dATACpeaks.V717I[dATACpeaks.V717I$annotation != "Promoter"]
#11,587 peaks
length(dATACpeaks.V717I.nonpromoter)


#Prepare promoter regions for APP-V717I 
#da_log2FC = Fold
V717I.dATAC.promoter.up <- dATACpeaks.V717I.promoter %>%
  dplyr::filter(da_log2FC > 0)
V717I.dATAC.promoter.down <- dATACpeaks.V717I.promoter %>%
  dplyr::filter(da_log2FC < 0)

#1,226 peaks
length(V717I.dATAC.promoter.up)

#3,839 peaks
length(V717I.dATAC.promoter.down)


#Convert granges object to a dataframe
V717I.dATAC.promoter.up.df <- as.data.frame(V717I.dATAC.promoter.up)
#Grab the seqnames, start, end 
V717I.dATAC.promoter.up.df <- V717I.dATAC.promoter.up.df[,c(1:3)]

#Convert granges object to a dataframe
V717I.dATAC.promoter.down.df <- as.data.frame(V717I.dATAC.promoter.down)
#Grab the seqnames, start, end
V717I.dATAC.promoter.down.df <- V717I.dATAC.promoter.down.df[,c(1:3)]

#############################################
# Run chipenrich on all promoter regions 
#############################################

#Use logistic regression model test function chipenrich in the chipenrich R package with the 
#locus definition nearest_tss for enrichment of promoter-located peaks and 1kb_outside for 
#enhancer-located peaks using the ENCODE-ChEA and ReMap TF gene-target databases and Gene 
#Onology (GO) - Biological Process and Hallmark ontology databases. 

#GOBP
#yes thats correct, the idea is that you take your filtered dATAC peaks that are promoter down, 
#and the create a dataframe that is a bed file format (first 3 columns of granges object)

#GOBP- Promoter Down Region
V717I.dATAC.promoter.down.GOBP <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.down.df,
                                                         genome = 'hg38',
                                                         genesets = "GOBP",
                                                         locusdef = "nearest_tss",
                                                         qc_plots = FALSE,
                                                         out_name = NULL,
                                                         n_cores=12)


V717I.dATAC.promoter.down.GOBP.results = V717I.dATAC.promoter.down.GOBP$results
print(V717I.dATAC.promoter.down.GOBP.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.down.GOBP.results, file = "~/ChIPenrich/V717I.dATAC.promoter.down.GOBP.chipenrich.results.csv")


#GOBP- Promoter Up Region
V717I.dATAC.promoter.up.GOBP <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.up.df,
                                                       genome = 'hg38',
                                                       genesets = "GOBP",
                                                       locusdef = "nearest_tss",
                                                       qc_plots = FALSE,
                                                       out_name = NULL,
                                                       n_cores=12)

V717I.dATAC.promoter.up.GOBP.results = V717I.dATAC.promoter.up.GOBP$results
print(V717I.dATAC.promoter.up.GOBP.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.up.GOBP.results, file = "~/ChIPenrich/V717I.dATAC.promoter.up.GOBP.chipenrich.results.csv")

#ENCODE-ChEA Consensus - Promoter Down Region
V717I.dATAC.promoter.down.ECC <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.down.df,
                                                        genome = 'hg38',
                                                        genesets = "~/gmtfiles/ECC_TFchipenrich.txt",
                                                        locusdef = "nearest_tss",
                                                        qc_plots = FALSE,
                                                        out_name = NULL,
                                                        n_cores=12)
V717I.dATAC.promoter.down.ECC.results = V717I.dATAC.promoter.down.ECC$results
print(V717I.dATAC.promoter.down.ECC.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.down.ECC.results, file = "~/ChIPenrich/V717I.dATAC.promoter.down.ECC.chipenrich.results.csv")

#ENCODE-ChEA Consensus - Promoter Up Region
V717I.dATAC.promoter.up.ECC <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.up.df,
                                                      genome = 'hg38',
                                                      genesets = "~/gmtfiles/ECC_TFchipenrich.txt",
                                                      locusdef = "nearest_tss",
                                                      qc_plots = FALSE,
                                                      out_name = NULL,
                                                      n_cores=12)
V717I.dATAC.promoter.up.ECC.results = V717I.dATAC.promoter.up.ECC$results
print(V717I.dATAC.promoter.up.ECC.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.up.ECC.results, file = "~/ChIPenrich/V717I.dATAC.promoter.up.ECC.chipenrich.results.csv")

#
#ReMap - Promoter Down Region
V717I.dATAC.promoter.down.ReMap <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.down.df,
                                                          genome = 'hg38',
                                                          genesets = "~/gmtfiles/ReMap_TFchipenrich.txt",
                                                          locusdef = "nearest_tss",
                                                          qc_plots = FALSE,
                                                          out_name = NULL,
                                                          n_cores=12)
V717I.dATAC.promoter.down.ReMap.results = V717I.dATAC.promoter.down.ReMap$results
print(V717I.dATAC.promoter.down.ReMap.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.down.ReMap.results, file = "~/ChIPenrich/V717I.dATAC.promoter.down.ReMap.chipenrich.results.csv")

#ReMap - Promoter Up Region
V717I.dATAC.promoter.up.ReMap <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.up.df,
                                                        genome = 'hg38',
                                                        genesets = "~/gmtfiles/ReMap_TFchipenrich.txt",
                                                        locusdef = "nearest_tss",
                                                        qc_plots = FALSE,
                                                        out_name = NULL,
                                                        n_cores=12)
V717I.dATAC.promoter.up.ReMap.results = V717I.dATAC.promoter.up.ReMap$results
print(V717I.dATAC.promoter.up.ReMap.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.up.ReMap.results, file = "~/ChIPenrich/V717I.dATAC.promoter.up.ReMap.chipenrich.results.csv")

#Hallmark - Promoter Down Region
V717I.dATAC.promoter.down.Hallmark <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.down.df,
                                                             genome = 'hg38',
                                                             genesets = "hallmark",
                                                             locusdef = "nearest_tss",
                                                             qc_plots = FALSE,
                                                             out_name = NULL,
                                                             n_cores=12)
V717I.dATAC.promoter.down.Hallmark.results = V717I.dATAC.promoter.down.Hallmark$results
print(V717I.dATAC.promoter.down.Hallmark.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.down.Hallmark.results, file = "~/ChIPenrich/V717I.dATAC.promoter.down.Hallmark.chipenrich.results.csv")

#Hallmark - Promoter Up Region
V717I.dATAC.promoter.up.Hallmark <- chipenrich::chipenrich(peaks = V717I.dATAC.promoter.up.df,
                                                           genome = 'hg38',
                                                           genesets = "hallmark",
                                                           locusdef = "nearest_tss",
                                                           qc_plots = FALSE,
                                                           out_name = NULL,
                                                           n_cores=12)
V717I.dATAC.promoter.up.Hallmark.results = V717I.dATAC.promoter.up.Hallmark$results
print(V717I.dATAC.promoter.up.Hallmark.results[1:50,2:5])
write.csv(V717I.dATAC.promoter.up.Hallmark.results, file = "~/ChIPenrich/V717I.dATAC.promoter.up.Hallmark.chipenrich.results.csv")


##################################################################################
# Running the ChIPenrich program using R (A79V vs. NDC)
##################################################################################

#28,772 peaks found
dATACpeaks.A79V <- anno.A79V.filtered

#Load the differential ATAC peaks for PSEN1-A79V

#28,772 peaks
length(dATACpeaks.A79V)

dATACpeaks.A79V.promoter <- dATACpeaks.A79V[dATACpeaks.A79V$annotation == "Promoter"]
#7,595 peaks
length(dATACpeaks.A79V.promoter)

dATACpeaks.A79V.nonpromoter <- dATACpeaks.A79V[dATACpeaks.A79V$annotation != "Promoter"]
#21,177 peaks
length(dATACpeaks.A79V.nonpromoter)


#Prepare promoter regions for PSEN1-A79V
#da_log2FC = Fold
A79V.dATAC.promoter.up <- dATACpeaks.A79V.promoter %>%
  dplyr::filter(da_log2FC > 0)
A79V.dATAC.promoter.down <- dATACpeaks.A79V.promoter %>%
  dplyr::filter(da_log2FC < 0)

#2,435 peaks
length(A79V.dATAC.promoter.up)

#5,160 peaks
length(A79V.dATAC.promoter.down)


#Convert granges object to a dataframe
A79V.dATAC.promoter.up.df <- as.data.frame(A79V.dATAC.promoter.up)
#Grab the seqnames, start, end 
A79V.dATAC.promoter.up.df <- A79V.dATAC.promoter.up.df[,c(1:3)]

#Convert granges object to a dataframe
A79V.dATAC.promoter.down.df <- as.data.frame(A79V.dATAC.promoter.down)
#Grab the seqnames, start, end
A79V.dATAC.promoter.down.df <- A79V.dATAC.promoter.down.df[,c(1:3)]

#############################################
# Run chipenrich on all promoter regions
#############################################

#Use logistic regression model test function chipenrich in the chipenrich R package with the 
#locus definition nearest_tss for enrichment of promoter-located peaks and 1kb_outside for 
#enhancer-located peaks using the ENCODE-ChEA and ReMap TF gene-target databases and Gene 
#Onology (GO) - Biological Process and Hallmark ontology databases. 

#GOBP
#yes thats correct, the idea is that you take your filtered dATAC peaks that are promoter down, 
#and the create a dataframe that is a bed file format (first 3 columns of granges object)


#GOBP- Promoter Down Region
A79V.dATAC.promoter.down.GOBP <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.down.df,
                                                        genome = 'hg38',
                                                        genesets = "GOBP",
                                                        locusdef = "nearest_tss",
                                                        qc_plots = FALSE,
                                                        out_name = NULL,
                                                        n_cores=12)


A79V.dATAC.promoter.down.GOBP.results = A79V.dATAC.promoter.down.GOBP$results
print(A79V.dATAC.promoter.down.GOBP.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.down.GOBP.results, file = "~/ChIPenrich/A79V.dATAC.promoter.down.GOBP.chipenrich.results.csv")


#GOBP- Promoter Up Region
A79V.dATAC.promoter.up.GOBP <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.up.df,
                                                      genome = 'hg38',
                                                      genesets = "GOBP",
                                                      locusdef = "nearest_tss",
                                                      qc_plots = FALSE,
                                                      out_name = NULL,
                                                      n_cores=12)

A79V.dATAC.promoter.up.GOBP.results = A79V.dATAC.promoter.up.GOBP$results
print(A79V.dATAC.promoter.up.GOBP.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.up.GOBP.results, file = "~/ChIPenrich/A79V.dATAC.promoter.up.GOBP.chipenrich.results.csv")

#ENCODE-ChEA Consensus - Promoter Down Region
A79V.dATAC.promoter.down.ECC <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.down.df,
                                                       genome = 'hg38',
                                                       genesets = "~/gmtfiles/ECC_TFchipenrich.txt",
                                                       locusdef = "nearest_tss",
                                                       qc_plots = FALSE,
                                                       out_name = NULL,
                                                       n_cores=12)
A79V.dATAC.promoter.down.ECC.results = A79V.dATAC.promoter.down.ECC$results
print(A79V.dATAC.promoter.down.ECC.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.down.ECC.results, file = "~/ChIPenrich/A79V.dATAC.promoter.down.ECC.chipenrich.results.csv")

#ENCODE-ChEA Consensus - Promoter Up Region
A79V.dATAC.promoter.up.ECC <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.up.df,
                                                     genome = 'hg38',
                                                     genesets = "~/gmtfiles/ECC_TFchipenrich.txt",
                                                     locusdef = "nearest_tss",
                                                     qc_plots = FALSE,
                                                     out_name = NULL,
                                                     n_cores=12)
A79V.dATAC.promoter.up.ECC.results = A79V.dATAC.promoter.up.ECC$results
print(A79V.dATAC.promoter.up.ECC.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.up.ECC.results, file = "~/ChIPenrich/A79V.dATAC.promoter.up.ECC.chipenrich.results.csv")

#
#ReMap - Promoter Down Region
A79V.dATAC.promoter.down.ReMap <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.down.df,
                                                         genome = 'hg38',
                                                         genesets = "~/gmtfiles/ReMap_TFchipenrich.txt",
                                                         locusdef = "nearest_tss",
                                                         qc_plots = FALSE,
                                                         out_name = NULL,
                                                         n_cores=12)
A79V.dATAC.promoter.down.ReMap.results = A79V.dATAC.promoter.down.ReMap$results
print(A79V.dATAC.promoter.down.ReMap.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.down.ReMap.results, file = "~/ChIPenrich/A79V.dATAC.promoter.down.ReMap.chipenrich.results.csv")

#ReMap - Promoter Up Region
A79V.dATAC.promoter.up.ReMap <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.up.df,
                                                       genome = 'hg38',
                                                       genesets = "~/gmtfiles/ReMap_TFchipenrich.txt",
                                                       locusdef = "nearest_tss",
                                                       qc_plots = FALSE,
                                                       out_name = NULL,
                                                       n_cores=12)
A79V.dATAC.promoter.up.ReMap.results = A79V.dATAC.promoter.up.ReMap$results
print(A79V.dATAC.promoter.up.ReMap.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.up.ReMap.results, file = "~/ChIPenrich/A79V.dATAC.promoter.up.ReMap.chipenrich.results.csv")

#Hallmark - Promoter Down Region
A79V.dATAC.promoter.down.Hallmark <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.down.df,
                                                            genome = 'hg38',
                                                            genesets = "hallmark",
                                                            locusdef = "nearest_tss",
                                                            qc_plots = FALSE,
                                                            out_name = NULL,
                                                            n_cores=12)
A79V.dATAC.promoter.down.Hallmark.results = A79V.dATAC.promoter.down.Hallmark$results
print(A79V.dATAC.promoter.down.Hallmark.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.down.Hallmark.results, file = "~/ChIPenrich/A79V.dATAC.promoter.down.Hallmark.chipenrich.results.csv")

#Hallmark - Promoter Up Region
A79V.dATAC.promoter.up.Hallmark <- chipenrich::chipenrich(peaks = A79V.dATAC.promoter.up.df,
                                                          genome = 'hg38',
                                                          genesets = "hallmark",
                                                          locusdef = "nearest_tss",
                                                          qc_plots = FALSE,
                                                          out_name = NULL,
                                                          n_cores=12)
A79V.dATAC.promoter.up.Hallmark.results = A79V.dATAC.promoter.up.Hallmark$results
print(A79V.dATAC.promoter.up.Hallmark.results[1:50,2:5])
write.csv(A79V.dATAC.promoter.up.Hallmark.results, file = "~/ChIPenrich/A79V.dATAC.promoter.up.Hallmark.chipenrich.results.csv")


####################################################################################
# Running the ChIPenrich program using R (N141I vs. NDC)
####################################################################################


#19,721 peaks found
dATACpeaks.N141I <- anno.N141I.filtered

#Load the differential ATAC peaks for PSEN2-N141I

#19,721 peaks
length(dATACpeaks.N141I)

dATACpeaks.N141I.promoter <- dATACpeaks.N141I[dATACpeaks.N141I$annotation == "Promoter"]
#5,543 peaks
length(dATACpeaks.N141I.promoter)

dATACpeaks.N141I.nonpromoter <- dATACpeaks.N141I[dATACpeaks.N141I$annotation != "Promoter"]
#14,781 peaks
length(dATACpeaks.N141I.nonpromoter)


#Prepare promoter regions for PSEN2-N141I
#da_log2FC = Fold
N141I.dATAC.promoter.up <- dATACpeaks.N141I.promoter %>%
  dplyr::filter(da_log2FC > 0)
N141I.dATAC.promoter.down <- dATACpeaks.N141I.promoter %>%
  dplyr::filter(da_log2FC < 0)

#1,547 peaks
length(N141I.dATAC.promoter.up)

#3,996 peaks
length(N141I.dATAC.promoter.down)


#Convert granges object to a dataframe
N141I.dATAC.promoter.up.df <- as.data.frame(N141I.dATAC.promoter.up)
#Grab the seqnames, start, end 
N141I.dATAC.promoter.up.df <- N141I.dATAC.promoter.up.df[,c(1:3)]

#Convert granges object to a dataframe
N141I.dATAC.promoter.down.df <- as.data.frame(N141I.dATAC.promoter.down)
#Grab the seqnames, start, end
N141I.dATAC.promoter.down.df <- N141I.dATAC.promoter.down.df[,c(1:3)]

#############################################
# Run chipenrich on all promoter regions
#############################################

#Use logistic regression model test function chipenrich in the chipenrich R package with the 
#locus definition nearest_tss for enrichment of promoter-located peaks and 1kb_outside for 
#enhancer-located peaks using the ENCODE-ChEA and ReMap TF gene-target databases and Gene 
#Onology (GO) - Biological Process and Hallmark ontology databases. 

#GOBP
#yes thats correct, the idea is that you take your filtered dATAC peaks that are promoter down, 
#and the create a dataframe that is a bed file format (first 3 columns of granges object)

#GOBP- Promoter Down Region
N141I.dATAC.promoter.down.GOBP <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.down.df,
                                                         genome = 'hg38',
                                                         genesets = "GOBP",
                                                         locusdef = "nearest_tss",
                                                         qc_plots = FALSE,
                                                         out_name = NULL,
                                                         n_cores=12)


N141I.dATAC.promoter.down.GOBP.results = N141I.dATAC.promoter.down.GOBP$results
print(N141I.dATAC.promoter.down.GOBP.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.down.GOBP.results, file = "~/ChIPenrich/N141I.dATAC.promoter.down.GOBP.chipenrich.results.csv")


#GOBP- Promoter Up Region
N141I.dATAC.promoter.up.GOBP <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.up.df,
                                                       genome = 'hg38',
                                                       genesets = "GOBP",
                                                       locusdef = "nearest_tss",
                                                       qc_plots = FALSE,
                                                       out_name = NULL,
                                                       n_cores=12)

N141I.dATAC.promoter.up.GOBP.results = N141I.dATAC.promoter.up.GOBP$results
print(N141I.dATAC.promoter.up.GOBP.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.up.GOBP.results, file = "~/ChIPenrich/N141I.dATAC.promoter.up.GOBP.chipenrich.results.csv")

#ENCODE-ChEA Consensus - Promoter Down Region
N141I.dATAC.promoter.down.ECC <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.down.df,
                                                        genome = 'hg38',
                                                        genesets = "~/gmtfiles/ECC_TFchipenrich.txt",
                                                        locusdef = "nearest_tss",
                                                        qc_plots = FALSE,
                                                        out_name = NULL,
                                                        n_cores=12)
N141I.dATAC.promoter.down.ECC.results = N141I.dATAC.promoter.down.ECC$results
print(N141I.dATAC.promoter.down.ECC.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.down.ECC.results, file = "~/ChIPenrich/N141I.dATAC.promoter.down.ECC.chipenrich.results.csv")

#ENCODE-ChEA Consensus - Promoter Up Region
N141I.dATAC.promoter.up.ECC <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.up.df,
                                                      genome = 'hg38',
                                                      genesets = "~/gmtfiles/ECC_TFchipenrich.txt",
                                                      locusdef = "nearest_tss",
                                                      qc_plots = FALSE,
                                                      out_name = NULL,
                                                      n_cores=12)
N141I.dATAC.promoter.up.ECC.results = N141I.dATAC.promoter.up.ECC$results
print(N141I.dATAC.promoter.up.ECC.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.up.ECC.results, file = "~/ChIPenrich/N141I.dATAC.promoter.up.ECC.chipenrich.results.csv")

#
#ReMap - Promoter Down Region
N141I.dATAC.promoter.down.ReMap <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.down.df,
                                                          genome = 'hg38',
                                                          genesets = "~/gmtfiles/ReMap_TFchipenrich.txt",
                                                          locusdef = "nearest_tss",
                                                          qc_plots = FALSE,
                                                          out_name = NULL,
                                                          n_cores=12)
N141I.dATAC.promoter.down.ReMap.results = N141I.dATAC.promoter.down.ReMap$results
print(N141I.dATAC.promoter.down.ReMap.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.down.ReMap.results, file = "~/ChIPenrich/N141I.dATAC.promoter.down.ReMap.chipenrich.results.csv")

#ReMap - Promoter Up Region
N141I.dATAC.promoter.up.ReMap <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.up.df,
                                                        genome = 'hg38',
                                                        genesets = "~/gmtfiles/ReMap_TFchipenrich.txt",
                                                        locusdef = "nearest_tss",
                                                        qc_plots = FALSE,
                                                        out_name = NULL,
                                                        n_cores=12)
N141I.dATAC.promoter.up.ReMap.results = N141I.dATAC.promoter.up.ReMap$results
print(N141I.dATAC.promoter.up.ReMap.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.up.ReMap.results, file = "~/ChIPenrich/N141I.dATAC.promoter.up.ReMap.chipenrich.results.csv")


#Hallmark - Promoter Down Region
N141I.dATAC.promoter.down.Hallmark <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.down.df,
                                                             genome = 'hg38',
                                                             genesets = "hallmark",
                                                             locusdef = "nearest_tss",
                                                             qc_plots = FALSE,
                                                             out_name = NULL,
                                                             n_cores=12)
N141I.dATAC.promoter.down.Hallmark.results = N141I.dATAC.promoter.down.Hallmark$results
print(N141I.dATAC.promoter.down.Hallmark.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.down.Hallmark.results, file = "~N141I.dATAC.promoter.down.Hallmark.chipenrich.results.csv")

#Hallmark - Promoter Up Region
N141I.dATAC.promoter.up.Hallmark <- chipenrich::chipenrich(peaks = N141I.dATAC.promoter.up.df,
                                                           genome = 'hg38',
                                                           genesets = "hallmark",
                                                           locusdef = "nearest_tss",
                                                           qc_plots = FALSE,
                                                           out_name = NULL,
                                                           n_cores=12)
N141I.dATAC.promoter.up.Hallmark.results = N141I.dATAC.promoter.up.Hallmark$results
print(N141I.dATAC.promoter.up.Hallmark.results[1:50,2:5])
write.csv(N141I.dATAC.promoter.up.Hallmark.results, file = "~/ChIPenrich/N141I.dATAC.promoter.up.Hallmark.chipenrich.results.csv")

#########################################
# Save the .RData File for chipenrich
#########################################
save.image('~/ATAC-seq/RData/ChIPEnrich_Env.RData')
