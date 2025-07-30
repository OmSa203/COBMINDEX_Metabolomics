#----------------Cleaning the environment
rm(list=ls())


library("ggplot2")
library("dplyr")
library("IRdisplay")
library("readxl")
library("plyr")

library(tidyverse)
library(factoextra)
library(KODAMA)
library(vegan)



#Function for saving. if the directory doesn't exist it will create new.
check_dir_and_save_csv <- function(Data_to_save, directory, filename) {
  # Construct full directory path
  full_dir <- directory
  
  # Create directory if it doesn't exist
  if (!dir.exists(full_dir)) {
    dir.create(full_dir, recursive = TRUE)
    message("Created directory: ", full_dir)
  }
  
  # Construct the full file path
  final_file_path <- file.path(full_dir, filename)
  
  # Write the CSV
  write.csv(Data_to_save, file = final_file_path, row.names = FALSE)
  return()
}

#The plot saving doesnt work- TO CHECK
Frequency_plots <- function(freqTable,plotTitle){
  FreqTable <- freqTable
  pTitle <- plotTitle
  final_plot <- ggplot(FreqTable, aes(x=Range_Bins, y=Log_Freq)) + 
    geom_bar(stat = "identity", position = "dodge", width=0.3) + 
    ggtitle(label = pTitle) +
    xlab("Range") + 
    ylab("(Log)Frequency") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  Cutoff_LOD <- round(min(blk_rem[blk_rem > 0]))
  message(paste0("The limit of detection (LOD) is: ",Cutoff_LOD)) 
  return(final_plot)
}
#===========================For Negative========================================

#shows working directory
getwd()
#set the R working directory
setwd("C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data")
#show files in directory
dir()
#Reading R Data Files+ naming the files

#Quant_Table output
FT <- read.csv("NEW_Quant_Table_Negative.csv")
head(FT)

#metadata
MD <- read.csv("NEGmetaData.csv")
head(MD)


#Selecting the columns with samples
New_FT <- select(FT, matches("mzML"))

#Removing "peak area" from column names (to match MD)
colnames(New_FT) <- gsub(".Peak.area","",colnames(New_FT))

#Changing names in metadata (to match New_FT)
MD$filename <- gsub('-','.',MD$filename)

#Finding blanks
Blank <- filter(MD, Sample_Type == "Blank")
Blank_data <- select(New_FT, c(Blank$filename))

#Finding samples
Sample <- filter(MD, Sample_Type == "Sample")
Sample_data <- select(New_FT, c(Sample$filename))

#Adding average to blanks and samples and adding row.ID
Blank_data$MeanB <-apply(Blank_data,1,mean)
Sample_data$MeanS <-apply(Sample_data,1,mean)
Blank_data <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z), RT = select(FT, row.retention.time), Blank_data)
Sample_data <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z), RT = select(FT, row.retention.time), Sample_data)

#Combining averages into new file
BSaverage <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z),RT = select(FT, row.retention.time), MeanB = select(Blank_data, MeanB), MeanS = select(Sample_data, MeanS))

#Blank/Sample ratio
BSaverage$ratio <- (BSaverage$MeanB+1)/(BSaverage$MeanS+1)

#if ratio is bigger than 2 then it's 0, if not then it's 1
BSaverage <- BSaverage %>% 
  mutate(bin = if_else(BSaverage$ratio > 0.3, 1, 0))

#Number of features after blank removal
FeaturesRemoved <- (nrow(BSaverage)-sum(BSaverage$bin))
FeaturesRemained <- sum(BSaverage$bin)

#Filtering all samples and blank removal
Sample_data <- subset (Sample_data, select = -MeanS)
rownames(Sample_data) <- Sample_data[,1]
Sample2 <- filter(BSaverage, bin == "1")
Sample2 <- data.frame(Sample2[,-1], row.names=Sample2[,1])
Sample_clean <- Sample_data %>% filter(rownames(Sample_data) %in% rownames(Sample2))

save_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Cleaned_data/"
check_dir_and_save_csv(Sample_clean,save_dir,"Data_Cleanup_Negative.csv")


#==============Imputation ================

blk_rem <- t(subset(Sample_clean, select = -c(2,3)))
colnames(blk_rem) <- blk_rem[1,]  # Set column names to first row
blk_rem <- blk_rem[-1,]  # Remove the first row
blk_rem <- as.data.frame(blk_rem)

#creating bins from -1 to 10^10 using sequence function seq()
bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

#cut function cuts the give table into its appropriate bins
scores_gapfilled <- cut(as.matrix(blk_rem),bins, labels = c('0','1',paste("1E",1:10,sep="")))

#transform function convert the tables into a column format: easy for visualization
FreqTable <- transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values
colnames(FreqTable)[1] <- 'Range_Bins' #changing the 1st colname to 'Range Bins'

## Frequency Plot
freq_plot <- Frequency_plots(FreqTable,"Frequency plot - Gap to Fill: Negative")
freq_plot

setwd(save_dir)
ggsave("Negative Frequency plot - Gap to Fill.png", freq_plot, width = 15, height = 10)



set.seed(141222) # by setting a seed, we generate the same set of random number all the time

imp <- blk_rem

for (i in 1:ncol(imp)) {
  imp[,i] <- ifelse(imp[,i] == 0, 
                    round(runif(nrow(imp), min = 1, max = Cutoff_LOD), 1),
                    imp[,i])
}

#creating bins from -1 to 10^10 using sequence function seq()
binsi <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

#cut function cuts the give table into its appropriate bins
scores_gapfilled <- cut(as.matrix(imp),binsi, labels = c('0','1',paste("1E",1:10,sep="")))

#transform function convert the tables into a column format: easy for visualization
FreqTable <- transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values
colnames(FreqTable)[1] <- 'Range_Bins' #changing the 1st colname to 'Range Bins'

## Frequency Plot
freq_plot_filled <- Frequency_plots(FreqTable,"Frequency plot - Gap Filled: Negative")
freq_plot_filled
ggsave("Negative Frequency plot - Gap Filled.png", freq_plot_filled, width = 15, height = 10)


sum(imp==0) # checking if there are any zeros in our imputed table 
dir()
imp2 <- as.data.frame(t(imp))
imp_Data <- data.frame(row.ID = select(Sample_clean, row.ID), 
                       mz = select(Sample_clean, row.m.z), 
                       RT = select(Sample_clean, row.retention.time)
                       ,imp2)
write.csv(blk_rem, 'Before_Imputed_QuantTable_Negative.csv', row.names =TRUE)
write.csv(imp_Data, 'Imputed_QuantTable_Negative.csv', row.names =FALSE)


#==============Normalization=============================
##Total Ion Current (TIC) or sample-centric normalization

norm_tic <- normalization(imp, 
                          method = "sum")$newXtrain
head(norm_tic,n=3)
dim(norm_tic)
print(paste('No.of NA values in Normalized data:', sum(is.na(norm_tic)== T)))

norm_tic <- as.data.frame(t(norm_tic))
Norma_Data <- data.frame(row.ID = select(Sample_clean, row.ID), 
                         mz = select(Sample_clean, row.m.z), 
                         RT = select(Sample_clean, row.retention.time)
                         ,norm_tic)

write.csv(Norma_Data, 'Imputed_TIC_Normalised_table_Negative.csv', row.names =FALSE)





#===========================For Positive========================================
#----------------Cleaning the environment
rm(list=ls())



#**CALL the check_dir_and_save_csv function again**
#shows working directory
getwd()
#set the R working directory
setwd("C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/Data")
#show files in directory
dir()
#Reading R Data Files+ naming the files

#Quant_Table output
FT <- read.csv("NEW_Quant_Table_Positive.csv")
head(FT)

#metadata
MD <- read.csv("POSmetaData.csv")
head(MD)


#Selecting the columns with samples
New_FT <- select(FT, matches("mzML"))

#Removing "peak area" from column names (to match MD)
colnames(New_FT) <- gsub(".Peak.area","",colnames(New_FT))

#Changing names in metadata (to match New_FT)
MD$filename <- gsub('-','.',MD$filename)

#Finding blanks
Blank <- filter(MD, Sample_Type == "Blank")
Blank_data <- select(New_FT, c(Blank$filename))

#Finding samples
Sample <- filter(MD, Sample_Type == "Sample")
Sample_data <- select(New_FT, c(Sample$filename))

#Adding average to blanks and samples and adding row.ID
Blank_data$MeanB <-apply(Blank_data,1,mean)
Sample_data$MeanS <-apply(Sample_data,1,mean)
Blank_data <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z), RT = select(FT, row.retention.time), Blank_data)
Sample_data <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z), RT = select(FT, row.retention.time), Sample_data)

#Combining averages into new file
BSaverage <- data.frame(row.ID = select(FT, row.ID), mz = select(FT, row.m.z),RT = select(FT, row.retention.time), MeanB = select(Blank_data, MeanB), MeanS = select(Sample_data, MeanS))

#Blank/Sample ratio
BSaverage$ratio <- (BSaverage$MeanB+1)/(BSaverage$MeanS+1)

#if ratio is bigger than 2 then it's 0, if not then it's 1
BSaverage <- BSaverage %>% 
  mutate(bin = if_else(BSaverage$ratio > 0.3, 1, 0))

#Number of features after blank removal
FeaturesRemoved <- (nrow(BSaverage)-sum(BSaverage$bin))
FeaturesRemained <- sum(BSaverage$bin)

#Filtering all samples and blank removal
Sample_data <- subset (Sample_data, select = -MeanS)
rownames(Sample_data) <- Sample_data[,1]
Sample2 <- filter(BSaverage, bin == "1")
Sample2 <- data.frame(Sample2[,-1], row.names=Sample2[,1])
Sample_clean <- Sample_data %>% filter(rownames(Sample_data) %in% rownames(Sample2))

save_dir <- "C:/Users/omers/Documents/GitHub/COBMINDEX_Metabolomics/results/Cleaned_data/"
check_dir_and_save_csv(Sample_clean,save_dir,"Data_Cleanup_Positive.csv")

#==============Imputation ================

blk_rem <- t(subset(Sample_clean, select = -c(2,3)))
colnames(blk_rem) <- blk_rem[1,]  # Set column names to first row
blk_rem <- blk_rem[-1,]  # Remove the first row
blk_rem <- as.data.frame(blk_rem)

#creating bins from -1 to 10^10 using sequence function seq()
bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

#cut function cuts the give table into its appropriate bins
scores_gapfilled <- cut(as.matrix(blk_rem),bins, labels = c('0','1',paste("1E",1:10,sep="")))

#transform function convert the tables into a column format: easy for visualization
FreqTable <- transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values
colnames(FreqTable)[1] <- 'Range_Bins' #changing the 1st colname to 'Range Bins'

## Frequency Plot
freq_plot <- Frequency_plots(FreqTable,"Frequency plot - Gap Filled: Positive")
freq_plot


set.seed(141222) # by setting a seed, we generate the same set of random number all the time

imp <- blk_rem

for (i in 1:ncol(imp)) {
  imp[,i] <- ifelse(imp[,i] == 0, 
                    round(runif(nrow(imp), min = 1, max = Cutoff_LOD), 1),
                    imp[,i])
}

#creating bins from -1 to 10^10 using sequence function seq()
binsi <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

#cut function cuts the give table into its appropriate bins
scores_gapfilled <- cut(as.matrix(imp),binsi, labels = c('0','1',paste("1E",1:10,sep="")))

#transform function convert the tables into a column format: easy for visualization
FreqTable <- transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values
colnames(FreqTable)[1] <- 'Range_Bins' #changing the 1st colname to 'Range Bins'

## Frequency Plot
freq_plot_filled <- Frequency_plots(FreqTable,"Imputed Frequency plot - Gap Filled: Positive")
freq_plot_filled

setwd(save_dir)
ggsave("Imputed Frequency plot - Gap Filled: Positive.png", freq_plot_filled, width = 15, height = 10)

sum(imp==0) # checking if there are any zeros in our imputed table 
dir()
imp2 <- as.data.frame(t(imp))
imp_Data <- data.frame(row.ID = select(Sample_clean, row.ID), 
                         mz = select(Sample_clean, row.m.z), 
                         RT = select(Sample_clean, row.retention.time)
                         ,imp2)

write.csv(blk_rem, 'Before_Imputed_QuantTable_Positive.csv', row.names =TRUE)
write.csv(imp_Data, 'Imputed_QuantTable_Positive.csv', row.names =FALSE)


#==================Normalization=============================


##Total Ion Current (TIC) or sample-centric normalization

norm_tic <- normalization(imp, 
                          method = "sum")$newXtrain
head(norm_tic,n=3)
dim(norm_tic)
print(paste('No.of NA values in Normalized data:', sum(is.na(norm_tic)== T)))

norm_tic <- as.data.frame(t(norm_tic))
Norma_Data <- data.frame(row.ID = select(Sample_clean, row.ID), 
                         mz = select(Sample_clean, row.m.z), 
                         RT = select(Sample_clean, row.retention.time)
                         ,norm_tic)

write.csv(Norma_Data, "Imputed_TIC_Normalised_table_Positive.csv", row.names =FALSE)
dir()

sessionInfo()

