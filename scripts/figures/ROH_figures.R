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

#adding group column
nroh_sroh$group <- sub("_.*", "", nroh_sroh$samples)

#Create a Scatter plot 
library(ggplot2)

ggplot(nroh_sroh, aes(x = TOTAL / 1000000, y = nTOTAL, color = group)) +  # Scale X by 1 million
  geom_point(size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +  # Add x = y line
  labs(title = "ROH",
       x = "TOTAL (in millions of bp)",  # Adjust X-axis label to reflect scale
       y = "nTOTAL",  # Y-axis stays the same
       color = "Group") +
  scale_color_manual(values = c("BBIP" = "#61bab8", "SBJO" = "#deef7f", "SJVE" = "#faa56f")) +  # Custom colors
  theme_light()




  



