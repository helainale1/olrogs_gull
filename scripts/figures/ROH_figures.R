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

#Create a Scatter plot 
library(ggplot2)

ggplot(nroh_sroh, aes(x = TOTAL / 1000000, y = nTOTAL, color = group)) +  # Scale X by 1 million
  geom_point(size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +  # Add x = y line
  labs(title = "ROH",
       x = "TOTAL (in millions)",  # Adjust X-axis label to reflect scale
       y = "nTOTAL (hundreds)",  # Y-axis stays the same
       color = "Group") +
  scale_color_manual(values = c("BBIP" = "#817fc7", "SBJO" = "#fe829b", "SJVE" = "#40e0d0")) +  # Custom colors
  theme_light()




  



