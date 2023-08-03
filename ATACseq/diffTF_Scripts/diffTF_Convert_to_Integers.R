--------------------------------------------------------------------------------------------------------
#Create raw counts table for input into diffTF 
  
#Create filtered raw counts data table
final_counts = as.data.frame(y$counts)
  
#Move the index columns into the first column called ENSEMBL
final_counts <- cbind(ENSEMBL = rownames(final_counts), final_counts)
View(final_counts)
  
#Source: https://stackoverflow.com/questions/6081439/changing-column-names-of-a-data-frame
#Rename names of colummns
colnames(final_counts) <- c("ENSEMBL","NDC_1", "NDC_2", "NDC_3", "V717I_1", "V717I_2", "V717I_3",
                              "A79V_1", "A79V_2", "A79V_3", "N141I_1", "N141I_2", "N141I_3")
  
final_counts = as.data.frame(final_counts)
  
#Convert counts to integers
final_counts_num = as.matrix(final_counts[2:13])
class(final_counts_num) <- "numeric"
  
#Source: https://stackoverflow.com/questions/36943953/fastest-way-to-coerce-matrix-to-integer-matrix-in-r
#Convert to integer
mode(final_counts_num) <- "integer"
  
#Change the rownames
rownames(final_counts_num) <- final_counts$ENSEMBL
  
#Move the index columns into the first column called ENSEMBL
final_counts_num <- cbind(ENSEMBL = rownames(final_counts_num), final_counts_num)
  
 #Create new .csv file
write.csv(final_counts_num, 
            file="/Users/phoebevaldes/Desktop/diffTF/Counts/RNA-seq-counts2.csv", 
            sep="\t", quote=F, row.names=F, col.names = T)