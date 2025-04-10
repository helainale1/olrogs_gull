# Load ggplot2 library
library(ggplot2)
roh_summary <- read.table("/Users/hdl5108/desktop/olrogs_gull/data/olrogs_roh_sum_0.5MB.txt")

#Rename rows
samples <- rownames(roh_summary)
roh_summary$samples <- samples
rownames(roh_summary) <- NULL

# Remove "L_atlanticus_" from the 'samples' column
roh_summary$samples <- gsub("^L_atlanticus_", "", roh_summary$samples)

# Remove the pop column 
roh_summary$pop <- NULL 

#create new data frame
nroh_sroh <- roh_summary[, c("samples", "nTOTAL", "TOTAL")]

