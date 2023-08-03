#ATAC-seq/GWAS integration Analysis pipeline 
#Valdes et. al 2023 Molecular Psychiatry Submission
#Load packages
library("plyranges")
library("granges")

###################################################################################################
#14. Perform integration of AD GWAS loci and differential accessible regions (DARs) at FDR < 0.05
#APP-V717I vs. NDC
###################################################################################################

#####################################################################
#14a. AD GWAS loci Table from compiled literature and database resources
#####################################################################

#Read the compiled AD GWAS loci study table
ADloci_Table = read.csv("~/GWAS_Loci_Analysis/AD_variants_chr_pos_Table.csv")

#Convert AD GWAS loci study table to GRanges object
ADloci_Table_Granges <- plyranges::as_granges(ADloci_Table)

####################################
#14b. APP-V717I vs. NDC in all DARs
####################################

#Find overlap between the two granges objects (ATACSeq.report.V717I and ADloci_Table_Granges)
AD_loci_APPvs.NDC_GRanges_Overlap <- plyranges::join_overlap_inner(ATACSeq.report.V717I, ADloci_Table_Granges)

#Get unique terms
AD_loci_APPvs.NDC_GRanges_Overlap <- unique(AD_loci_APPvs.NDC_GRanges_Overlap)

#Convert to dataframe
AD_loci_APPvs.NDC_GRanges_Overlap_df <- as.data.frame(AD_loci_APPvs.NDC_GRanges_Overlap)

#Write to .csv file
write.csv(AD_loci_APPvs.NDC_GRanges_Overlap_df, file="~/GWAS_Loci_Analysis/APPvsNDC_AD_variants_DARs.csv", sep="\t", quote=F, row.names=T)

#####################################################
#14c. APP-V717I vs. NDC in promoter-associated DARs
#####################################################

#Find overlap between the two granges objects (dATACpeaks.V717I.promoter and ADloci_Table_Granges)
AD_loci_APPvs.NDC_GRanges_Overlap.promoter <- plyranges::join_overlap_inner(dATACpeaks.V717I.promoter, ADloci_Table_Granges)

#Get unique terms
AD_loci_APPvs.NDC_GRanges_Overlap.promoter <- unique(AD_loci_APPvs.NDC_GRanges_Overlap.promoter)

#Convert to dataframe
AD_loci_APPvs.NDC_GRanges_Overlap_promoter_df <- as.data.frame(AD_loci_APPvs.NDC_GRanges_Overlap.promoter)

#Write to .csv file
write.csv(AD_loci_APPvs.NDC_GRanges_Overlap_promoter_df, file="~/GWAS_Loci_Analysis/APPvsNDC_AD_variants_DARs_promoter.csv", sep="\t", quote=F, row.names=T)

#################################################################################
#14d. APP-V717I vs. NDC in promoter-associated DARs with increased accessibility
#################################################################################

#Find overlap between the two granges objects (V717I.dATAC.promoter.up and ADloci_Table_Granges)
AD_loci_APPvs.NDC_GRanges_Overlap.promoter.up <- plyranges::join_overlap_inner(V717I.dATAC.promoter.up, ADloci_Table_Granges)

#Get unique terms
AD_loci_APPvs.NDC_GRanges_Overlap.promoter.up <- unique(AD_loci_APPvs.NDC_GRanges_Overlap.promoter.up)

#Convert to dataframe
AD_loci_APPvs.NDC_GRanges_Overlap_promoter_up_df <- as.data.frame(AD_loci_APPvs.NDC_GRanges_Overlap.promoter.up)

#Write to .csv file
write.csv(AD_loci_APPvs.NDC_GRanges_Overlap_promoter_up_df, file="~/GWAS_Loci_Analysis/APPvsNDC_AD_variants_DARs_promoter_up.csv", sep="\t", quote=F, row.names=T)

#################################################################################
#14e. APP-V717I vs. NDC in promoter-associated DARs with decreased accessibility
#################################################################################

#Find overlap between the two granges objects (V717I.dATAC.promoter.down and ADloci_Table_Granges)
AD_loci_APPvs.NDC_GRanges_Overlap.promoter.down <- plyranges::join_overlap_inner(V717I.dATAC.promoter.down, ADloci_Table_Granges)

#Get unique terms
AD_loci_APPvs.NDC_GRanges_Overlap.promoter.down <- unique(AD_loci_APPvs.NDC_GRanges_Overlap.promoter.down)

#Convert to dataframe
AD_loci_APPvs.NDC_GRanges_Overlap_promoter_down_df <- as.data.frame(AD_loci_APPvs.NDC_GRanges_Overlap.promoter.down)

#Write to .csv file
write.csv(AD_loci_APPvs.NDC_GRanges_Overlap_promoter_down_df, file="~/GWAS_Loci_Analysis/APPvsNDC_AD_variants_DARs_promoter_down.csv", sep="\t", quote=F, row.names=T)


###############################################################
#14f. APP-V717I vs. NDC in PEREGRINE enhancer-associated DARs
###############################################################

#Find overlap between the two granges objects (dATACpeaks.V717I.enhancer.PEREGRINE and ADloci_Table_Granges)
AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer <- plyranges::join_overlap_inner(dATACpeaks.V717I.enhancer.PEREGRINE, ADloci_Table_Granges)

#Get unique terms
AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer <- unique(AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer)

#Convert to dataframe
AD_loci_APPvs.NDC_GRanges_Overlap_PEREGRINE.enhancer_df <- as.data.frame(AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer)

#Write to .csv file
write.csv(AD_loci_APPvs.NDC_GRanges_Overlap_PEREGRINE.enhancer_df, file="~/GWAS_Loci_Analysis/APPvsNDC_AD_variants_DARs_PEREGRINE.enhancer.csv", sep="\t", quote=F, row.names=T)

###########################################################################################
#14g. APP-V717I vs. NDC in PEREGRINE enhancer-associated DARs with increased accessibility
###########################################################################################

#Find overlap between the two granges objects (V717I.dATAC.enhancer.PEREGRINE.up and ADloci_Table_Granges)
AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up <- plyranges::join_overlap_inner(V717I.dATAC.enhancer.PEREGRINE.up, ADloci_Table_Granges)

#Get unique terms
AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up <- unique(AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up)

#Convert to dataframe
AD_loci_APPvs.NDC_GRanges_Overlap_PEREGRINE.enhancer_up_df <- as.data.frame(AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up)

#Write to .csv file
write.csv(AD_loci_APPvs.NDC_GRanges_Overlap_PEREGRINE.enhancer_up_df, file="~/GWAS_Loci_Analysis/APPvsNDC_AD_variants_DARs_PEREGRINE.enhancer_up.csv", sep="\t", quote=F, row.names=T)

###########################################################################################
#14h. APP-V717I vs. NDC in PEREGRINE enhancer-associated DARs with decreased accessibility
###########################################################################################

#Find overlap between the two granges objects (V717I.dATAC.enhancer.PEREGRINE.down and ADloci_Table_Granges)
AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down <- plyranges::join_overlap_inner(V717I.dATAC.enhancer.PEREGRINE.down, ADloci_Table_Granges)

#Get unique terms
AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down <- unique(AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down)

#Convert to dataframe
AD_loci_APPvs.NDC_GRanges_Overlap_PEREGRINE.enhancer_down_df <- as.data.frame(AD_loci_APPvs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down)

#Write to .csv file
write.csv(AD_loci_APPvs.NDC_GRanges_Overlap_PEREGRINE.enhancer_down_df, file="~/GWAS_Loci_Analysis/APPvsNDC_AD_variants_DARs_PEREGRINE.enhancer_down.csv", sep="\t", quote=F, row.names=T)

###################################################################################################
#15. Perform integration of AD GWAS loci and differential accessible regions (DARs) at FDR < 0.05
#PSEN1-A79V vs. NDC
###################################################################################################

#########################
#15a. PSEN1-A79V vs. NDC
#########################

#Find overlap between the two granges objects (ATACSeq.report.A79V and ADloci_Table_Granges)
AD_loci_PSEN1vs.NDC_GRanges_Overlap <- plyranges::join_overlap_inner(ATACSeq.report.A79V, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN1vs.NDC_GRanges_Overlap <- unique(AD_loci_PSEN1vs.NDC_GRanges_Overlap)

#Convert to dataframe
AD_loci_PSEN1vs.NDC_GRanges_Overlap_df <- as.data.frame(AD_loci_PSEN1vs.NDC_GRanges_Overlap)

#Write to .csv file
write.csv(AD_loci_PSEN1vs.NDC_GRanges_Overlap_df, file="~/GWAS_Loci_Analysis/PSEN1vsNDC_variants/230608_PSEN1vs.NDC_AD_variants_DARs.csv", sep="\t", quote=F, row.names=T)

#####################################################
#15b. PSEN1-A79V vs. NDC in promoter-associated DARs
#####################################################

#Find overlap between the two granges objects (dATACpeaks.A79V.promoter and ADloci_Table_Granges)
AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter <- plyranges::join_overlap_inner(dATACpeaks.A79V.promoter, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter <- unique(AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter)

#Convert to dataframe
AD_loci_PSEN1vs.NDC_GRanges_Overlap_promoter_df <- as.data.frame(AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter)

#Write to .csv file
write.csv(AD_loci_PSEN1vs.NDC_GRanges_Overlap_promoter_df, file="~/GWAS_Loci_Analysis/PSEN1vsNDC_AD_variants_DARs_promoter.csv", sep="\t", quote=F, row.names=T)

##################################################################################
#15c. PSEN1-A79V vs. NDC in promoter-associated DARs with increased accessibility
##################################################################################

#Find overlap between the two granges objects (A79V.dATAC.promoter.up and ADloci_Table_Granges)
AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.up <- plyranges::join_overlap_inner(A79V.dATAC.promoter.up, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.up <- unique(AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.up)

#Convert to dataframe
AD_loci_PSEN1vs.NDC_GRanges_Overlap_promoter_up_df <- as.data.frame(AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.up)

#Write to .csv file
write.csv(AD_loci_PSEN1vs.NDC_GRanges_Overlap_promoter_up_df, file="~/GWAS_Loci_Analysis/PSEN1vsNDC_AD_variants_DARs_promoter_up.csv", sep="\t", quote=F, row.names=T)

##################################################################################
#15d. PSEN1-A79V vs. NDC in promoter-associated DARs with decreased accessibility
##################################################################################

#Find overlap between the two granges objects (A79V.dATAC.promoter.down and ADloci_Table_Granges)
AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.down <- plyranges::join_overlap_inner(A79V.dATAC.promoter.down, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.down <- unique(AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.down)

#Convert to dataframe
AD_loci_PSEN1vs.NDC_GRanges_Overlap_promoter_down_df <- as.data.frame(AD_loci_PSEN1vs.NDC_GRanges_Overlap.promoter.down)

#Write to .csv file
write.csv(AD_loci_PSEN1vs.NDC_GRanges_Overlap_promoter_down_df, file="~/GWAS_Loci_Analysis/PSEN1vsNDC_AD_variants_DARs_promoter_down.csv", sep="\t", quote=F, row.names=T)


###############################################################
#15e. PSEN1-A79V vs. NDC in PEREGRINE enhancer-associated DARs
###############################################################

#Find overlap between the two granges objects (dATACpeaks.A79V.enhancer.PEREGRINE and ADloci_Table_Granges)
AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer <- plyranges::join_overlap_inner(dATACpeaks.A79V.enhancer.PEREGRINE, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer <- unique(AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer)

#Convert to dataframe
AD_loci_PSEN1vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_df <- as.data.frame(AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer)

#Write to .csv file
write.csv(AD_loci_PSEN1vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_df, file="~/GWAS_Loci_Analysis/PSEN1vsNDC_AD_variants_DARs_PEREGRINE.enhancer.csv", sep="\t", quote=F, row.names=T)

############################################################################################
#15f. PSEN1-A79V vs. NDC in PEREGRINE enhancer-associated DARs with increased accessibility
############################################################################################

#Find overlap between the two granges objects (A79V.dATAC.enhancer.PEREGRINE.up and ADloci_Table_Granges)
AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up <- plyranges::join_overlap_inner(A79V.dATAC.enhancer.PEREGRINE.up, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up <- unique(AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up)

#Convert to dataframe
AD_loci_PSEN1vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_up_df <- as.data.frame(AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up)

#Write to .csv file
write.csv(AD_loci_PSEN1vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_up_df, file="~/GWAS_Loci_Analysis/PSEN1vsNDC_AD_variants_DARs_PEREGRINE.enhancer_up.csv", sep="\t", quote=F, row.names=T)

############################################################################################
#15g. PSEN1-A79V vs. NDC in PEREGRINE enhancer-associated DARs with decreased accessibility
############################################################################################

#Find overlap between the two granges objects (A79V.dATAC.enhancer.PEREGRINE.down and ADloci_Table_Granges)
AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down <- plyranges::join_overlap_inner(A79V.dATAC.enhancer.PEREGRINE.down, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down <- unique(AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down)

#Convert to dataframe
AD_loci_PSEN1vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_down_df <- as.data.frame(AD_loci_PSEN1vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down)

#Write to .csv file
write.csv(AD_loci_PSEN1vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_down_df, file="~/GWAS_Loci_Analysis/PSEN1vsNDC_AD_variants_DARs_PEREGRINE.enhancer_down.csv", sep="\t", quote=F, row.names=T)



###################################################################################################
#16. Perform integration of AD GWAS loci and differential accessible regions (DARs) at FDR < 0.05
#PSEN2-N141I vs. NDC
###################################################################################################

##########################
#16a. PSEN2-N141I vs. NDC
##########################

#Find overlap between the two granges objects (ATACSeq.report.N141I and ADloci_Table_Granges)
AD_loci_PSEN2vs.NDC_GRanges_Overlap <- plyranges::join_overlap_inner(ATACSeq.report.N141I, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN2vs.NDC_GRanges_Overlap <- unique(AD_loci_PSEN2vs.NDC_GRanges_Overlap)

#Convert to dataframe
AD_loci_PSEN2vs.NDC_GRanges_Overlap_df <- as.data.frame(AD_loci_PSEN2vs.NDC_GRanges_Overlap)

#Write to .csv file
write.csv(AD_loci_PSEN2vs.NDC_GRanges_Overlap_df, file="~/GWAS_Loci_Analysis/PSEN2vsNDC_AD_variants_DARs.csv", sep="\t", quote=F, row.names=T)


######################################################
#16b. PSEN2-N141I vs. NDC in promoter-associated DARs
######################################################

#Find overlap between the two granges objects (dATACpeaks.A79V.promoter and ADloci_Table_Granges)
AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter <- plyranges::join_overlap_inner(dATACpeaks.N141I.promoter, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter <- unique(AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter)

#Convert to dataframe
AD_loci_PSEN2vs.NDC_GRanges_Overlap_promoter_df <- as.data.frame(AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter)

#Write to .csv file
write.csv(AD_loci_PSEN2vs.NDC_GRanges_Overlap_promoter_df, file="~/GWAS_Loci_Analysis/PSEN2vsNDC_AD_variants_DARs_promoter.csv", sep="\t", quote=F, row.names=T)

###################################################################################
#16c. PSEN2-N141I vs. NDC in promoter-associated DARs with increased accessibility
###################################################################################

#Find overlap between the two granges objects (dATACpeaks.V717I.promoter and ADloci_Table_Granges)
AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.up <- plyranges::join_overlap_inner(N141I.dATAC.promoter.up, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.up <- unique(AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.up)

#Convert to dataframe
AD_loci_PSEN2vs.NDC_GRanges_Overlap_promoter_up_df <- as.data.frame(AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.up)

#Write to .csv file
write.csv(AD_loci_PSEN2vs.NDC_GRanges_Overlap_promoter_up_df, file="~/GWAS_Loci_Analysis/PSEN2vsNDC_AD_variants_DARs_promoter_up.csv", sep="\t", quote=F, row.names=T)


###################################################################################
#16d. PSEN2-N141I vs. NDC in promoter-associated DARs with decreased accessibility
###################################################################################

#Find overlap between the two granges objects (dATACpeaks.V717I.promoter and ADloci_Table_Granges)
AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.down <- plyranges::join_overlap_inner(N141I.dATAC.promoter.down, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.down <- unique(AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.down)

#Convert to dataframe
AD_loci_PSEN2vs.NDC_GRanges_Overlap_promoter_down_df <- as.data.frame(AD_loci_PSEN2vs.NDC_GRanges_Overlap.promoter.down)

#Write to .csv file
write.csv(AD_loci_PSEN2vs.NDC_GRanges_Overlap_promoter_down_df, file="~/GWAS_Loci_Analysis/PSEN2vsNDC_AD_variants_DARs_promoter_down.csv", sep="\t", quote=F, row.names=T)

################################################################
#16e. PSEN2-N141I vs. NDC in PEREGRINE enhancer-associated DARs
################################################################

#Find overlap between the two granges objects (dATACpeaks.N141I.enhancer.PEREGRINE and ADloci_Table_Granges)
AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer <- plyranges::join_overlap_inner(dATACpeaks.N141I.enhancer.PEREGRINE, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer <- unique(AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer)

#Convert to dataframe
AD_loci_PSEN2vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_df <- as.data.frame(AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer)

#Write to .csv file
write.csv(AD_loci_PSEN2vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_df, file="~/GWAS_Loci_Analysis/PSEN2vsNDC_AD_variants_DARs_PEREGRINE.enhancer.csv", sep="\t", quote=F, row.names=T)

#############################################################################################
#16f. PSEN2-N141I vs. NDC in PEREGRINE enhancer-associated DARs with increased accessibility
#############################################################################################

#Find overlap between the two granges objects (N141I.dATAC.enhancer.PEREGRINE.up and ADloci_Table_Granges)
AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up <- plyranges::join_overlap_inner(N141I.dATAC.enhancer.PEREGRINE.up, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up <- unique(AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up)

#Convert to dataframe
AD_loci_PSEN2vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_up_df <- as.data.frame(AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.up)

#Write to .csv file
write.csv(AD_loci_PSEN2vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_up_df, file="~/GWAS_Loci_Analysis/PSEN2vsNDC_AD_variants_DARs_PEREGRINE.enhancer_up.csv", sep="\t", quote=F, row.names=T)


##############################################################################################
#16g. PSEN2-N141I vs. NDC in PEREGRINE enhancer-associated DARs with decreased accessibility
##############################################################################################

#Find overlap between the two granges objects (N141I.dATAC.enhancer.PEREGRINE.down and ADloci_Table_Granges)
AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down <- plyranges::join_overlap_inner(N141I.dATAC.enhancer.PEREGRINE.down, ADloci_Table_Granges)

#Get unique terms
AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down <- unique(AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down)

#Convert to dataframe
AD_loci_PSEN2vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_down_df <- as.data.frame(AD_loci_PSEN2vs.NDC_GRanges_Overlap.PEREGRINE.enhancer.down)

#Write to .csv file
write.csv(AD_loci_PSEN2vs.NDC_GRanges_Overlap_PEREGRINE.enhancer_down_df, file="~/GWAS_Loci_Analysis/PSEN2vsNDC_AD_variants_DARs_PEREGRINE.enhancer_down.csv", sep="\t", quote=F, row.names=T)

###########################################################################################################################################################################################

#############################################
#17. Look at common variants among mutations
#############################################

#Get APP vs. NDC variants
APP_variants <- read.csv("~/GWAS_Loci_Analysis/APPvs.NDC_variants/APPvs.NDC_overlap_GWAS_variants.csv")
APP_variants_list <- list(APP_variants = as.character(unique(APP_variants$Variant)))


#Get PSEN1 vs. NDC variants  
PSEN1_variants <- read.csv("~/GWAS_Loci_Analysis/PSEN1vs.NDC_variants/PSEN1vs.NDC_overlap_GWAS_variants.csv")
PSEN1_variants_list <- list(PSEN1_variants = as.character(unique(PSEN1_variants$Variant)))


#Get PSEN2 vs. NDC variants  
PSEN2_variants <- read.csv("~/GWAS_Loci_Analysis/PSEN2vs.NDC_variants/PSEN2vs.NDC_overlap_GWAS_variants.csv")
PSEN2_variants_list <- list(PSEN2_variants = as.character(unique(PSEN2_variants$Variant)))


#Get common variants by variant ID
APP_PSEN1_common_variants <- merge(APP_variants, PSEN1_variants, by = "Variant")
PSEN1_PSEN2_common_variants <- merge(PSEN1_variants, PSEN2_variants, by = "Variant")
APP_PSEN2_common_variants <- merge(APP_variants, PSEN2_variants, by = "Variant")
APP_PSEN1_PSEN2_common_variants <- merge(APP_PSEN1_common_variants, PSEN2_variants, by = "Variant")


#Write .csv files
write.csv(APP_PSEN1_common_variants, file="~/GWAS_Loci_Analysis/Common_AD_variants/230621_Common_AD_variants_DARs_APP_PSEN1.csv", sep="\t", quote=F, row.names=T)
write.csv(PSEN1_PSEN2_common_variants, file="~/GWAS_Loci_Analysis/Common_AD_variants/230621_Common_AD_variants_DARs_PSEN1_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN2_common_variants, file="~/GWAS_Loci_Analysis/Common_AD_variants/230621_Common_AD_variants_DARs_APP_PSEN2.csv", sep="\t", quote=F, row.names=T)
write.csv(APP_PSEN1_PSEN2_common_variants, file="~/GWAS_Loci_Analysis/Common_AD_variants/230621_Common_AD_variants_DARs_allMutations.csv", sep="\t", quote=F, row.names=T)

#################################################################
#18. Create Venn Diagram for AD variants in all DARs (FDR < 0.05)
#################################################################

#Plot the nVenn Diagram 
myNV <- plotVenn(list(APP_variants_list, PSEN1_variants_list, PSEN2_variants_list), sNames=c("APP", "PSEN1", "PSEN2"),showPlot = T,nCycles = 10000)
showSVG(myNV, opacity=0.3,outFile = "~/GWAS_Loci_Analysis/Venn_Diagrams/All3_PSEN1_2_APP_Variants_nVenn_v1_NEW.svg", setColors = c("#fe6100", "#dc267f", "#785ef0"))


########################################################################
#19. Save the .RData file for ATAC-seq DiffBind with GWAS loci analysis
########################################################################
save.image('~/RData/071023_ATAC_Diffbind_VennDiagrams_Pipeline.RData')

#========================================================================================================