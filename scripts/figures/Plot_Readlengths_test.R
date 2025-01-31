test_readlengths <- read.table("/Users/hdl5108/desktop/olrogs_gull/data/test_readlength.txt")
colnames(test_readlengths) <- "read_length"
# Load ggplot2 library
library(ggplot2)

#Figure out binwidth
range(test_readlengths)
# Create a histogram of a numeric column in your data
ggplot(data = test_readlengths, aes(x = read_length)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "Sample 1",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))
