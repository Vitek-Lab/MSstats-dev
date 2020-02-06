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
LFQBench_data_set_SN_report_corrected_sel_normalization <- read_delim("MP-LFC-TTOF-SN-Report.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
# make column Run
data <- LFQBench_data_set_SN_report_corrected_sel_normalization %>% 
  unite(Run, R.Condition, R.Replicate, sep=".", remove = F) 

head(data)

info <- unique(data[, c('Run', 'R.Condition', 'Measurement', 'R.Replicate')])
info

write.csv(info, 'MP-LFC-TTOF-annotation.csv')


# check whether input has all the requried columns
.check.input(data)

#################### Step 2: Data preprocessing #######################
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


# constant normalization based on internal standard proteins
standards <- c("P68133", "Q71U36", "P60709", "P04350", "P62937", "P04406",
               "P18206", "P27797", "P63261", "P62979", "P60866", "Q14103",
               "E9PAV3", "Q15233", "P0CG48", "P63173", "P62913", "Q53S24",
               "Q9GZZ7", "P18124", "P60953", "O15372", "P62280", "Q9UNX3", "P61077")

# extract all the house-keeping proteins
proteins <- unique(data$IDPicker.InferenceId)
standards.all <- NULL
for(i in 1:length(standards)){
  standards.all <- c(standards.all, proteins[grepl(standards[i], proteins)])
}

# calculate the normalization factor
norm <- data %>% 
  filter(IDPicker.InferenceId %in% standards.all)  %>%
  select(Run, FG.MS1PeakArea, FG.MS2PeakArea) %>%
  group_by(Run) %>% 
  summarise(med.MS1 = median(FG.MS1PeakArea, na.rm = TRUE), med.MS2 = median(FG.MS2PeakArea, na.rm = TRUE)) %>%
  mutate(diff.MS1 = med.MS1 - median(med.MS1), diff.MS2 = med.MS2 - median(med.MS2))

# constant normalization
data <- data %>%
  left_join(norm) %>%
  mutate(FG.NormalizedMS1PeakArea = FG.MS1PeakArea - diff.MS1, 
         FG.NormalizedMS2PeakArea = FG.MS2PeakArea - diff.MS2) %>%
  select(-med.MS1, -diff.MS1, -med.MS2, -diff.MS2)

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

#################### Step 2 : remove precursors with missing values #######################
perc <- 1.0 # percentage of the missing values to remove
runs <- unique(filtered_data_long$RUN)
# remove the precursor which have any missing value in MS1 or MS2 data
filtered_data_long <- filtered_data_long[ , `:=`( COUNT = .N ) , by = .(PROTEIN, LABEL, FEATURE)][COUNT >= length(runs)*perc]
filtered_data_long[ ,COUNT:=NULL]

nrow(filtered_data_long)

################### Step 3: Run testing on MS1, MS2 and both #######################
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

