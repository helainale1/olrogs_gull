# Load the reshape2 package (if not installed yet)
library(reshape2)

# Adding population column 
roh_summary_clean$population <- sub("_.*", "", roh_summary_clean$samples)

# Melt everything into long format
roh_long <- melt(roh_summary_clean, id.vars = "samples")

# Split into "segments"
segments <- roh_long[grep("^n", roh_long$variable), ]
segments$class <- gsub("^n", "", segments$variable)
colnames(segments)[3] <- "segments"
segments <- segments[, c("samples", "class", "segments")]

# Merge only the segments dataset
roh_summary_long <- segments

# Adding population column to the merged data
roh_summary_long$population <- sub("_.*", "", roh_summary_long$samples)

# Reordering columns
roh_summary_long <- roh_summary_long[, c("population", "samples", "class", "segments")]

# Replace 'class' columns with each data frame 
roh_summary_long$class <- replace(
  roh_summary_long$class,                   # Target the 'class' column
  roh_summary_long$class == "A",            # Condition to find "A" values
  ".5-1MB"                                  # Replacement value
)

roh_summary_long$class <- replace(
  roh_summary_long$class,                   # Target the 'class' column
  roh_summary_long$class == "B",            # Condition to find "B" values
  "1-2MB"                                  # Replacement value
)

roh_summary_long$class <- replace(
  roh_summary_long$class,                   # Target the 'class' column
  roh_summary_long$class == "C",            # Condition to find "C" values
  "2-3MB"                                  # Replacement value
)

roh_summary_long$class <- replace(
  roh_summary_long$class,                   # Target the 'class' column
  roh_summary_long$class == "D",            # Condition to find "D" values
  "3-4MB"                                  # Replacement value
)

roh_summary_long$class <- replace(
  roh_summary_long$class,                   # Target the 'class' column
  roh_summary_long$class == "E",            # Condition to find "E" values
  "4-5MB"                                  # Replacement value
)

roh_summary_long$class <- replace(
  roh_summary_long$class,                   # Target the 'class' column
  roh_summary_long$class == "F",            # Condition to find "F" values
  ">5MB"                                  # Replacement value
)

roh_summary_long$segments <- as.numeric(roh_summary_long$segments)
roh_summary_long$class <- factor(roh_summary_long$class,
                                 levels = c(".5-1MB", "1-2MB", "2-3MB", "3-4MB", "4-5MB", ">5MB"))

# Check the result
head(roh_summary_long)
roh_summary_long$segments

#Creating the plot
# Load ggplot2 package
library(ggplot2)

ggplot(roh_summary_long, aes(x = class, y = segments, group = samples, color = population)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("BBIP" = "#61bab8", "SBJO" = "#deef7f", "SJVE" = "#faa56f")) +
  scale_y_continuous(breaks = seq(0,160, by=20), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") + 
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color="black", fill=NA),
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.margin = margin(0, 0, 0, 0)) + 
  labs(title = "Segments per Sample by Class", x = "Length of ROH", y = "Number of Segments") 

