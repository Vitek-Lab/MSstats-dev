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
PSSS2_Ting_classic_DIA_noMissingValues_v1.1 <- read_delim("Spike-in-biol-var-OT-SN-Report.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# make column Run
data <- PSSS2_Ting_classic_DIA_noMissingValues_v1.1 %>% 
  unite(Run, R.Condition, R.Replicate, sep=".", remove = F) 

head(data)

info <- unique(data[, c('Run', 'R.Condition', 'Measurement', 'R.Replicate')])
info

write.csv(info, 'Spike-in-biol-var-OT-annotation.csv')

# check whether input has all the requried columns
.check.input(data)

#################### Step 2: Data preprocessing #######################
# Remove the interfered blood proteins
rm.proteins <- c("P07724", "Q921I1", "ALBU_HUMAN_UPS")
data <- data %>% 
  filter(!IDPicker.InferenceId %in% rm.proteins)

# log transformation
data$FG.MS2PeakArea <- log2(data$FG.MS2PeakArea)
data$FG.NormalizedMS2PeakArea <- log2(data$FG.NormalizedMS2PeakArea)
data$FG.MS1PeakArea <- log2(data$FG.MS1PeakArea)
data$FG.NormalizedMS1PeakArea <- log2(data$FG.NormalizedMS1PeakArea)

# select the requried columns
data <- data %>% select(Run, R.Condition, R.Replicate, 
                        PG.Qvalue, EG.PrecursorId, EG.Qvalue, 
                        FG.MS2PeakArea, FG.NormalizedMS2PeakArea, 
                        FG.MS1PeakArea, FG.NormalizedMS1PeakArea, 
                        IDPicker.InferenceId)

# filter out the peptide ions whose MS1 signal <= 100 or qvalue > 0.01
filtered_data <- data %>%filter(EG.Qvalue <= 0.01 & 
                                  FG.NormalizedMS1PeakArea >= log2(100))
filtered_data <- as.data.table(filtered_data)

#################### Step 3: change data to long format #######################
filtered_data_long <- melt(filtered_data, 
                           id.vars = c("Run", "R.Condition", "R.Replicate", "IDPicker.InferenceId", "EG.PrecursorId"),
                           measure.vars = c("FG.NormalizedMS1PeakArea", "FG.NormalizedMS2PeakArea"), 
                           variable.name = "DataType", value.name = "PeakArea")
# rename the columns
setnames(filtered_data_long, 
         c("Run", "R.Condition", "R.Replicate", "IDPicker.InferenceId", "EG.PrecursorId", "DataType", "PeakArea"), 
         c("RUN", "GROUP", "SUBJECT", "PROTEIN", "FEATURE", "LABEL", "ABUNDANCE"))

nrow(filtered_data_long)

#################### Step 4: remove precursors with missing values #######################
perc <- 1.0 # percentage of the missing values to remove
runs <- unique(filtered_data_long$RUN)
# remove the precursor which have any missing value in MS1 or MS2 data
filtered_data_long <- filtered_data_long[ , `:=`( COUNT = .N ) , by = .(PROTEIN, LABEL, FEATURE)][COUNT >= length(runs)*perc]
filtered_data_long[ ,COUNT:=NULL]

nrow(filtered_data_long)

################### Step 5: Run testing on MS1, MS2 and both #######################
# fit linear model on the MS1 data
MS1.statistic <- protein.linear(data = filtered_data_long, MS = "FG.NormalizedMS1PeakArea")
save(MS1.statistic, file = "protein.linear.MS1.rda")
write.table(MS1.statistic , "NN-DIA protein.linear.MS1.txt", sep="\t", row.names = FALSE, quote = FALSE)
rm(MS1.statistic)

# fit linear model on the MS2 data
MS2.statistic <- protein.linear(data = filtered_data_long, MS = "FG.NormalizedMS2PeakArea")
save(MS2.statistic, file = "protein.linear.MS2.rda")
write.table(MS2.statistic , "NN-DIA protein.linear.MS2.txt", sep="\t", row.names = FALSE, quote = FALSE)
rm(MS2.statistic)

# fit linear model on the combined MS1 and MS2 data
MS1.MS2.statistic <- protein.linear(data = filtered_data_long, MS = c("FG.NormalizedMS1PeakArea", "FG.NormalizedMS2PeakArea"))
save(MS1.MS2.statistic, file = "protein.linear.MS1.MS2.rda")
write.table(MS1.MS2.statistic , "NN-DIA protein.linear.MS1.MS2.txt", sep="\t", row.names = FALSE, quote = FALSE)
rm(MS1.MS2.statistic)

