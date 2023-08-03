#ISMARA Report Extraction 
#Valdes et. al 2023 Molecular Psychiatry Submission
#Load packages

#Compute directional z-scores for V717I vs. NDC from ISMARA output
setwd("~/ISMARA_Results/V717IvsNDC.avg/averaged_report/")

#Function that extracts the ISMARA report for TF Regulon, Pearson score/p-value, z-value and activity difference
#for V717I vs. NDC
ismara_r_V717I <- function(y){
   active_matrix <- read.table("~/ISMARA_Results/V717IvsNDC.avg/averaged_report/active_matrices.txt", col.names = c("TF","Z value"),header = FALSE, stringsAsFactors = FALSE)
   filePath <- paste('~/ISMARA_Results/V717Ivs.NDCavg/averaged_report/correlation',active_matrix$TF,'correlation.html',sep='/')
   
      activity_matrix <- read.table("~/ISMARA_Results/V717IvsNDC.avg/averaged_report/activity_table.txt",header = TRUE, stringsAsFactors = FALSE)
      colnames(activity_matrix) <- gsub("\\.","-", colnames(activity_matrix))  #To make it uniform across all files

       file_sapply <- function(x) {  
           #Extract Motif's name
              unavailable <- strsplit(x,'/')
              unavailable <- unlist(unavailable)
              tff <- unavailable[8]
              
               #Extract Z-value of motif	 
               z_value <- active_matrix$Z.value[active_matrix$TF==tff[1]]
               
                  #Extract activities of the motifs 
                  activity_scores <- activity_matrix[,tff[1]]
                  
                   #Extract Pearson Co-efficient and the corresponding P-Value of various TFs that have the motif
                   if(file.exists(x)){
                       html_stuff <- readLines(x)
                       pearson_coef <- pearson(html_stuff)
                      }  else{
                         pearson_coef <- 0
                         return(0)
                       }
                 
                   pearson_coef$z_value <- z_value
                   pearson_coef$activity_diff <- activity_scores[2] - activity_scores[1]
                   
                   return(pearson_coef)
                 }
       
         pearson <- function(y){
               #Scrape the HTML page and isolate the line containing the necessary details - TF, Pearson Score, P-Value
               l <- length(y)
               l2 <- c(3:(l-1))
               t <- lapply(y[l2], data12)
               TFs <- data.frame(matrix(unlist(t), nrow = length(l2), byrow = TRUE), stringsAsFactors = FALSE)
               colnames(TFs) <- c("TFs", "Pearson_score", "Pearson_pvalue")
               return(TFs)
             }
         
           data12 <- function(y){
               #Extract the TF name, Pearson Score and P-Value
                 extract <- strsplit(y,'<td>')
                 extract <- unlist(extract)
                 TF <- gsub('(.*)</td>$','\\1',extract)[2]
                 Pearson <- as.numeric(gsub('(.*)</td>$','\\1',extract)[4])
                 Pearson_pvalue <- as.numeric(gsub('(.*)</td>$','\\1',extract)[5])
                 return(data.frame(TF, Pearson, Pearson_pvalue, stringsAsFactors = FALSE))
               }
           
             dat <- lapply(filePath, file_sapply)
             
               ismara_result <- do.call(rbind, dat)  #Merge sub-lists to form a data frame
               ismara_result <- ismara_result[ismara_result$TFs!=0,]
               head(ismara_result)
               names(ismara_result) <- c("Regulon",names(ismara_result)[c(2:length(names(ismara_result)))])
               fwrite(as.data.frame(ismara_result), file = file.path(paste0("~/ISMARA_Result","/V717IvsNDC"),paste0("ISMARA_",y,".csv")))
               return(ismara_result)
}

#Name the ISMARA .csv output as V717I_NDC
ismara_files <- as.list(c("V717I_NDC"))
ismara_results <- lapply(ismara_files, function(x){ismara_r_V717I(x)})

#Compute directional z-scores for A79V vs. NDC 
setwd("~/ISMARA_Results/A79VvsNDC.avg/averaged_report/")

#Function that extracts the ISMARA report for TF Regulon, Pearson score/p-value, z-value and activity difference
#for A79V vs. NDC
ismara_r_A79V <- function(y){
   active_matrix <- read.table("~/ISMARA_Results/A79VvsNDC.avg/averaged_report/active_matrices.txt", col.names = c("TF","Z value"),header = FALSE, stringsAsFactors = FALSE)
   filePath <- paste('~/ISMARA_Results/A79Vvs.NDC.avg/averaged_report/correlation',active_matrix$TF,'correlation.html',sep='/')
   
   activity_matrix <- read.table("~/ISMARA_Results/A79VvsNDC.avg/averaged_report/activity_table.txt",header = TRUE, stringsAsFactors = FALSE)
   colnames(activity_matrix) <- gsub("\\.","-", colnames(activity_matrix))  #To make it uniform across all files
   
   file_sapply <- function(x) {  
      #Extract Motif's name
      unavailable <- strsplit(x,'/')
      unavailable <- unlist(unavailable)
      tff <- unavailable[8]
      
      #Extract Z-value of motif	 
      z_value <- active_matrix$Z.value[active_matrix$TF==tff[1]]
      
      #Extract activities of the motifs 
      activity_scores <- activity_matrix[,tff[1]]
      
      #Extract Pearson Co-efficient and the corresponding P-Value of various TFs that have the motif
      if(file.exists(x)){
         html_stuff <- readLines(x)
         pearson_coef <- pearson(html_stuff)
      }  else{
         pearson_coef <- 0
         return(0)
      }
      
      pearson_coef$z_value <- z_value
      pearson_coef$activity_diff <- activity_scores[2] - activity_scores[1]
      
      return(pearson_coef)
   }
   
   pearson <- function(y){
      #Scrape the HTML page and isolate the line containing the necessary details - TF, Pearson Score, P-Value
      l <- length(y)
      l2 <- c(3:(l-1))
      t <- lapply(y[l2], data12)
      TFs <- data.frame(matrix(unlist(t), nrow = length(l2), byrow = TRUE), stringsAsFactors = FALSE)
      colnames(TFs) <- c("TFs", "Pearson_score", "Pearson_pvalue")
      return(TFs)
   }
   
   data12 <- function(y){
      #Extract the TF name, Pearson Score and P-Value
      extract <- strsplit(y,'<td>')
      extract <- unlist(extract)
      TF <- gsub('(.*)</td>$','\\1',extract)[2]
      Pearson <- as.numeric(gsub('(.*)</td>$','\\1',extract)[4])
      Pearson_pvalue <- as.numeric(gsub('(.*)</td>$','\\1',extract)[5])
      return(data.frame(TF, Pearson, Pearson_pvalue, stringsAsFactors = FALSE))
   }
   
   dat <- lapply(filePath, file_sapply)
   
   ismara_result <- do.call(rbind, dat)  #Merge sub-lists to form a data frame
   ismara_result <- ismara_result[ismara_result$TFs!=0,]
   head(ismara_result)
   names(ismara_result) <- c("Regulon",names(ismara_result)[c(2:length(names(ismara_result)))])
   fwrite(as.data.frame(ismara_result), file = file.path(paste0("~/ISMARA_Result","/A79VvsNDC"),paste0("ISMARA_",y,".csv")))
   return(ismara_result)
}

#Name the ISMARA .csv output as A79V_NDC
ismara_files <- as.list(c("A79V_NDC"))
ismara_results <- lapply(ismara_files, function(x){ismara_r_A79V(x)})


#Compute directional z-scores for N141I vs. NDC
setwd("~/ISMARA_Results/N141Ivs.NDC.avg/averaged_report/")

#Function that extracts the ISMARA report for TF Regulon, Pearson score/p-value, z-value and activity difference
#for N141I vs. NDC
ismara_r_N141I <- function(y){
   active_matrix <- read.table("~/ISMARA_Results/N141IvsNDC.avg/averaged_report/active_matrices.txt", col.names = c("TF","Z value"),header = FALSE, stringsAsFactors = FALSE)
   filePath <- paste('~/ISMARA_Results/N141IvsNDC.avg/averaged_report/correlation',active_matrix$TF,'correlation.html',sep='/')
   
   activity_matrix <- read.table("~/ISMARA_Results/N141IvsNDC.avg/averaged_report/activity_table.txt",header = TRUE, stringsAsFactors = FALSE)
   colnames(activity_matrix) <- gsub("\\.","-", colnames(activity_matrix))  #To make it uniform across all files
   
   file_sapply <- function(x) {  
      #Extract Motif's name
      unavailable <- strsplit(x,'/')
      unavailable <- unlist(unavailable)
      tff <- unavailable[8]
      
      #Extract Z-value of motif	 
      z_value <- active_matrix$Z.value[active_matrix$TF==tff[1]]
      
      #Extract activities of the motifs 
      activity_scores <- activity_matrix[,tff[1]]
      
      #Extract Pearson Co-efficient and the corresponding P-Value of various TFs that have the motif
      if(file.exists(x)){
         html_stuff <- readLines(x)
         pearson_coef <- pearson(html_stuff)
      }  else{
         pearson_coef <- 0
         return(0)
      }
      
      pearson_coef$z_value <- z_value
      pearson_coef$activity_diff <- activity_scores[2] - activity_scores[1]
      
      return(pearson_coef)
   }
   
   pearson <- function(y){
      #Scrape the HTML page and isolate the line containing the necessary details - TF, Pearson Score, P-Value
      l <- length(y)
      l2 <- c(3:(l-1))
      t <- lapply(y[l2], data12)
      TFs <- data.frame(matrix(unlist(t), nrow = length(l2), byrow = TRUE), stringsAsFactors = FALSE)
      colnames(TFs) <- c("TFs", "Pearson_score", "Pearson_pvalue")
      return(TFs)
   }
   
   data12 <- function(y){
      #Extract the TF name, Pearson Score and P-Value
      extract <- strsplit(y,'<td>')
      extract <- unlist(extract)
      TF <- gsub('(.*)</td>$','\\1',extract)[2]
      Pearson <- as.numeric(gsub('(.*)</td>$','\\1',extract)[4])
      Pearson_pvalue <- as.numeric(gsub('(.*)</td>$','\\1',extract)[5])
      return(data.frame(TF, Pearson, Pearson_pvalue, stringsAsFactors = FALSE))
   }
   
   dat <- lapply(filePath, file_sapply)
   
   ismara_result <- do.call(rbind, dat)  #Merge sub-lists to form a data frame
   ismara_result <- ismara_result[ismara_result$TFs!=0,]
   head(ismara_result)
   names(ismara_result) <- c("Regulon",names(ismara_result)[c(2:length(names(ismara_result)))])
   fwrite(as.data.frame(ismara_result), file = file.path(paste0("~/ISMARA_Result","/N141IvsNDC"),paste0("ISMARA_",y,".csv")))
   return(ismara_result)
}

#Find common TF regulons for all mutations from ISMARA results and then match by TF regulon

#Name the ISMARA .csv output as N141I_NDC
ismara_files <- as.list(c("N141I_NDC"))
ismara_results <- lapply(ismara_files, function(x){ismara_r_N141I(x)})

#Load packages
library("data.table")
library("gplots")
library("RColorBrewer")
library("dplyr")
library("tidyr")
library("ggplot2")

# Read in csv files
ISMARA_APP = read.csv('~/ISMARA_Results/ISMARA_V717I_NDC_dir_z-scores.csv')
ISMARA_PSEN1 = read.csv('~/ISMARA_Results/ISMARA_A79V_NDC_dir_z-scores.csv')
ISMARA_PSEN2 = read.csv('~/ISMARA_Results/ISMARA_N141I_NDC_dir_z-scores.csv')

#Convert each object into a data frame
ISMARA_APP = as.data.frame(ISMARA_APP)
ISMARA_PSEN1 = as.data.frame(ISMARA_PSEN1)
ISMARA_PSEN2 = as.data.frame(ISMARA_PSEN2)

#Filter for any TF regulons that meet the threshold of Pearson p-value < 0.01
ISMARA_APP.filtered <- ISMARA_APP %>%
  dplyr::filter(Pearson_pvalue < 0.1) 

ISMARA_PSEN1.filtered <- ISMARA_PSEN1 %>%
  dplyr::filter(Pearson_pvalue < 0.1) 

ISMARA_PSEN2.filtered <- ISMARA_PSEN2 %>%
  dplyr::filter(Pearson_pvalue < 0.1) 

write.csv(ISMARA_APP.filtered, file="~/ISMARA_Results/ISMARA_APP_filtered_Pearson-p-val.csv", row.names = FALSE)
write.csv(ISMARA_PSEN1.filtered, file="~/ISMARA_Results/ISMARA_PSEN1_filtered_Pearson-p-val.csv", row.names = FALSE)
write.csv(ISMARA_PSEN2.filtered, file="~/ISMARA_Results/ISMARA_PSEN2_filtered_Pearson-p-val.csv", row.names = FALSE)

############################################################################################
#Create filtered directional z-scores based on Pearson p-value <0.1
############################################################################################

#RNA-seq Figures 

#Read in .csv files
ISMARA.filtered_summary <- read.csv('~/ISMARA_Results/ISMARA_Summary_filtered_Pearson-p-val_dir_z-scores.csv')

#Convert to data frame
ISMARA.filtered_summary.df <- data.frame(ISMARA.filtered_summary)

#Convert variables into factors
ISMARA.filtered_summary.df$Contrast = as.factor(ISMARA.filtered_summary.df$Contrast)
new_colors_ISMARA.filtered_summary = c("APP-V717I vs. NDC", "PSEN1-A79V vs. NDC", "PSEN2-N141I vs. NDC")[ISMARA.filtered_summary.df$Contrast]

#Keep the ordering of the y axis, TF Regulon the 
#same as the data frame, ISMARA.filtered_summary
#based on the factor, Function (updated 4/19/21)

#Help Source: https://stackoverflow.com/questions/22340053/ggplot-order-of-factors-with-duplicate-levels

#Lock in factor level order by dataframe
ISMARA.filtered_summary.df$Regulon <- factor(ISMARA.filtered_summary.df$Regulon, levels = unique(ISMARA.filtered_summary.df$Regulon))


####################
#Merge TF Motifs
####################

#Regulon = Regulon
#Directional_Z-score.x = APP, Directional_Z-score.y = PSEN1 and Directional_Z_score = PSEN2
merged_TF_motifs = merge(merge(ISMARA_APP, ISMARA_PSEN1, by="Regulon"), ISMARA_PSEN2, by = "Regulon")

#All TF's have Pearson p-value < 0.10
merged_TF_motifs2 = merge(merge(ISMARA_APP.filtered, ISMARA_PSEN1.filtered, by="Regulon"), ISMARA_PSEN2.filtered, by = "Regulon")

#Create new data frame with merged directional z-scores 
#Source for help: https://stackoverflow.com/questions/50074987/create-new-dataframe-using-columns-from-multiple-data-frames
dir_z_scores_df = data.frame(merged_TF_motifs$Regulon, merged_TF_motifs$Directional_Z.score.x, merged_TF_motifs$Directional_Z.score.y, merged_TF_motifs$Directional.z.value)
dir_z_scores_df2 = data.frame(merged_TF_motifs2$Regulon, merged_TF_motifs2$Directional_Z.score.x, merged_TF_motifs2$Directional_Z.score.y, merged_TF_motifs2$Directional.z.value)


#Transpose the TF motif matrix
dir_z_scores_df <- as.data.frame(t(as.matrix(dir_z_scores_df)))
View(dir_z_scores_df)

dir_z_scores_df2 <- as.data.frame(t(as.matrix(dir_z_scores_df2)))
View(dir_z_scores_df2)

#Rename the columns
(setattr(dir_z_scores_df, "row.names", c("TF_Regulon", "APP-V717I", "PSEN1-A79V", "PSEN2-N141I")))
View(dir_z_scores_df)

(setattr(dir_z_scores_df2, "row.names", c("TF_Regulon", "APP-V717I", "PSEN1-A79V", "PSEN2-N141I")))
View(dir_z_scores_df2)

#Retranspose the TF motif matrix again
dir_z_scores_df <- t(dir_z_scores_df)
View(dir_z_scores_df)

#Retranspose the TF motif matrix again
dir_z_scores_df2 <- t(dir_z_scores_df2)
View(dir_z_scores_df2)

#Convert matrix into a data frame
dir_z_scores_df <- as.data.frame(dir_z_scores_df)
dir_z_scores_df2 <- as.data.frame(dir_z_scores_df2)


#Write merged TF motif results to .csv file
write.csv(dir_z_scores_df, file="~/ISMARA_Results/Merge_Files/merged_TF_regulons_dir_z_scores.csv", row.names = FALSE)  
write.csv(dir_z_scores_df2, file="~/ISMARA_Results/Merge_Files/merged_TF_regulons_dir_z_scores2.csv", row.names = FALSE)  

##############################
#Sort the merged TF regulons
##############################

#Read in the csv's after sorting by decreasing directional z-scores using Excel
dir_z_scores_df_sorted <- read.csv("~/ISMARA_Results/Merge_Files/merged_TF_regulons_dir_z_scores_sorted.csv")
dir_z_scores_df_sorted2 <- read.csv("~/ISMARA_Results/Merge_Files/merged_TF_regulons_dir_z_scores_sorted2.csv")

#For all the merged TF regulons
df_sorted_dir_z_scores_top <- dir_z_scores_df_sorted[1:40,] #top 40 regulons
df_sorted_dir_z_scores_top2 <- dir_z_scores_df_sorted[622:661,] #bottom 40 regulons

#For all the merged TF regulons that meet the Pearson score < 0.1
df_sorted_dir_z_scores_top_Pearson_p_val <- dir_z_scores_df_sorted2[1:40,]

head(df_sorted_dir_z_scores_top)
head(df_sorted_dir_z_scores_top2)
head(df_sorted_dir_z_scores_top_Pearson_p_val)

###################################################
#Extract the merged TF regulon as a numeric matrix
###################################################

#Extract just the numeric data into a matrix with named rows by gene
rownames(df_sorted_dir_z_scores_top) <- df_sorted_dir_z_scores_top$TF_Regulon
ISMARA_sorted_dir_z_scores_top <- as.matrix(df_sorted_dir_z_scores_top[2:4])

rownames(df_sorted_dir_z_scores_top2) <- df_sorted_dir_z_scores_top2$TF_Regulon
ISMARA_sorted_dir_z_scores_top2 <- as.matrix(df_sorted_dir_z_scores_top2[2:4])

rownames(df_sorted_dir_z_scores_top_Pearson_p_val) <- df_sorted_dir_z_scores_top_Pearson_p_val$TF_Regulon
ISMARA_sorted_dir_z_scores_Pearson_p_val <- as.matrix(df_sorted_dir_z_scores_top_Pearson_p_val[2:4])

#Remove the quotes from the matrix
class(ISMARA_sorted_dir_z_scores_top) <- "numeric"
class(ISMARA_sorted_dir_z_scores_top2) <- "numeric"
class(ISMARA_sorted_dir_z_scores_Pearson_p_val) <- "numeric"


