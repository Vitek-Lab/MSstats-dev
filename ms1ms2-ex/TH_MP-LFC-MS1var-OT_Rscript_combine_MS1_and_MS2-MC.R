library(readr)
library(tidyr)
library(dplyr)
library(data.table)
library(lme4)
library(nlme)

## read in the functions in LinearModelFunctions.R
source('LinearModelFunctions.R')

#################### Step 1: Read in the MS1 and MS2 data from Spectronaut #######################
# read in spectronaut peptide level data from both MS1 and MS2
data <- read_delim("MP-LFC-MS1var-OT-SN-Report.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(data)
unique(data$Measurement)

# extract the Run, R.Condition, R.Replicate information
data$Run <- data$Measurement
data$Run <- gsub("C_D170818_S441-", "", data$Run)
data$Run <- gsub("E_D170331_S209-", "", data$Run)
data$Run <- gsub("_T0", "", data$Run)
data <- data %>% separate(Run, c("R.Condition", "R.Replicate"), sep="_MHRM_", remove = FALSE)
data$Run <- gsub("_MHRM", "", data$Run)

info <- unique(data[, c('Run', 'R.Condition', 'Measurement', 'R.Replicate')])
info

write.csv(info, 'MP-LFC-MS1var-OT-annotation.csv')


#################### Step 2: Data preprocessing #######################
# log2 transformation
data$FG.MS2PeakArea <- log2(data$FG.MS2PeakArea)
data$FG.NormalizedMS2PeakArea <- as.numeric(data$FG.NormalizedMS2PeakArea)
data$FG.NormalizedMS2PeakArea <- log2(data$FG.NormalizedMS2PeakArea)
data$FG.MS1PeakArea <- log2(data$FG.MS1PeakArea)
data$FG.NormalizedMS1PeakArea <- as.numeric(data$FG.NormalizedMS1PeakArea)
data$FG.NormalizedMS1PeakArea <- log2(data$FG.NormalizedMS1PeakArea)

# select the requried columns
data <- data %>% 
  filter(!R.Condition %in% c("S-1-240min", "S-2-240min")) %>% 
  select(Run, R.Condition, R.Replicate, PG.Qvalue, 
         PG.ProteinAccessions, SpikeIn,
         PG.BGSOrganism, EG.PrecursorId, EG.Qvalue, 
         FG.MS2PeakArea, FG.NormalizedMS2PeakArea, 
         FG.MS1PeakArea, FG.NormalizedMS1PeakArea)

# constant normalization based on internal standard proteins
standards <- c("P68133", "Q71U36", "P60709", "P04350", "P62937", "P04406",
               "P18206", "P27797", "P63261", "P62979", "P60866", "Q14103",
               "E9PAV3", "Q15233", "P0CG48", "P63173", "P62913", "Q53S24",
               "Q9GZZ7", "P18124", "P60953", "O15372", "P62280", "Q9UNX3", "P61077")

# extract all the house-keeping proteins
proteins <- unique(data$PG.ProteinAccessions)
standards.all <- NULL
for(i in 1:length(standards)){
  standards.all <- c(standards.all, proteins[grepl(standards[i], proteins)])
}

# calculate the normalization factor for each dataset
norm <- data %>% 
  filter(PG.ProteinAccessions %in% standards.all)  %>%
  select(Run, FG.MS1PeakArea, FG.MS2PeakArea) %>%
  group_by(Run) %>% 
  summarise(med.MS1 = median(FG.MS1PeakArea, na.rm = TRUE), med.MS2 = median(FG.MS2PeakArea, na.rm = TRUE)) 

dt <- c("30k", "60k", "120k", "240k") # four datasets
norm2 <- list()
for(i in 1:4){
  temp <- norm %>% 
    filter(grepl(dt[i], Run))  %>%
    mutate(diff.MS1 = med.MS1 - median(med.MS1), diff.MS2 = med.MS2 - median(med.MS2))
  temp <- as.data.table(temp)
  norm2[[i]] <- temp 
}
norm2 <- do.call("rbind", norm2)

# constant normalization
data <- data %>%
  left_join(norm2) %>%
  mutate(FG.NormalizedMS1PeakArea = FG.MS1PeakArea - diff.MS1, 
       FG.NormalizedMS2PeakArea = FG.MS2PeakArea - diff.MS2) %>%
  select(-med.MS1, -diff.MS1, -med.MS2, -diff.MS2)

# filter out the peptide ions whose MS1 signal <= 100 or qvalue > 0.01
filtered_data <- data %>%filter(EG.Qvalue <= 0.01 & 
                                  FG.NormalizedMS1PeakArea >= log2(100))
filtered_data <- as.data.table(filtered_data)

#################### Step 3: statistical inference per dataset #######################
dt <- c("30k", "60k", "120k", "240k")
all.data <- filtered_data
for(i in 1:4){
  sub_data <- as.data.table(all.data %>% filter(grepl(dt[i], R.Condition)))
  
  sub_data_long <- melt(sub_data, 
                        id.vars = c("Run", "R.Condition", "R.Replicate", "PG.ProteinAccessions", "EG.PrecursorId"),
                        measure.vars = c("FG.NormalizedMS1PeakArea", "FG.NormalizedMS2PeakArea"), 
                        variable.name = "DataType", value.name = "PeakArea")
  
  # rename the columns
  setnames(sub_data_long, 
           c("Run", "R.Condition", "R.Replicate", "PG.ProteinAccessions", "EG.PrecursorId", "DataType", "PeakArea"), 
           c("RUN", "GROUP", "SUBJECT", "PROTEIN", "FEATURE", "LABEL", "ABUNDANCE"))
  
  nrow(sub_data_long)
  
  #################### Step 2 : remove precursors with missing values #######################
  perc <- 1.0 # percentage of the missing values to remove
  runs <- unique(sub_data_long$RUN)
  # remove the precursor which have any missing value in MS1 or MS2 data
  sub_data_long <- sub_data_long[ , `:=`( COUNT = .N ) , by = .(PROTEIN, LABEL, FEATURE)][COUNT >= length(runs)*perc]
  sub_data_long[ ,COUNT:=NULL]
  
  nrow(sub_data_long)
  
  ################### Step 3: Run testing on MS1, MS2 and both #######################
  # fit linear model on the MS1 data
  MS1.statistic <- protein.linear(data = sub_data_long, MS = "FG.NormalizedMS1PeakArea")
  
  save(MS1.statistic, file = paste(dt[i], ".protein.linear.MS1.rda", sep=""))
  write.table(MS1.statistic , file = paste("NN-DIA ", dt[i], ".protein.linear.MS1.txt", sep=""),
              sep="\t", row.names = FALSE, quote = FALSE)
  rm(MS1.statistic)
  
  # fit linear model on the MS2 data
  MS2.statistic <- protein.linear(data = sub_data_long, MS = "FG.NormalizedMS2PeakArea")
  save(MS2.statistic, file = paste(dt[i], ".protein.linear.MS2.rda", sep=""))
  write.table(MS2.statistic , file = paste("NN-DIA ", dt[i], ".protein.linear.MS2.txt", sep=""),
              sep="\t", row.names = FALSE, quote = FALSE)
  rm(MS2.statistic)
  
  # fit linear model on the combined MS1 and MS2 data
  MS1.MS2.statistic <- protein.linear(data = sub_data_long, MS = c("FG.NormalizedMS1PeakArea", "FG.NormalizedMS2PeakArea"))
  save(MS1.MS2.statistic, file = paste(dt[i], ".protein.linear.MS1.MS2.rda", sep=""))
  write.table(MS1.MS2.statistic , file = paste("NN-DIA ", dt[i], ".protein.linear.MS1.MS2.txt", sep=""),
              sep="\t", row.names = FALSE, quote = FALSE)
  rm(MS1.MS2.statistic)
}