# Plot Readlength Distributions

## Load R Module and Install Packages
You only have to install packages once
```bash
module use /storage/icds/RISE/sw8/modules/r 
module load r/4.2.1-gcc-8.5.0

Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('ggplot2', lib='/storage/group/dut374/default/bin/.R')"

Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('scales', lib='/storage/group/dut374/default/bin/.R')"

Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('ggpubr', lib='/storage/group/dut374/default/bin/.R')"
```

## Create R Script
```R
scripts_folder="/storage/group/dut374/default/helaina/scripts"
nano $scripts_folder/plots/plot_readlengths.R

#!/usr/bin/env Rscript

# CALL PACKAGES
library(ggplot2, lib.loc='/storage/group/dut374/default/bin/.R')
library(ggpubr, lib.loc='/storage/group/dut374/default/bin/.R')
library(scales, lib.loc='/storage/group/dut374/default/bin/.R')

#IMPORT DATA
BBIP_2_S1_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_BBIP_2_S1_readlength.txt")
BBIP_3_S2_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_BBIP_3_S2_readlength.txt")
BBIP_4_S3_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_BBIP_4_S3_readlength.txt")
BBIP_5_S4_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_BBIP_5_S4_readlength.txt")
BBIP_6_S5_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_BBIP_6_S5_readlength.txt")
SBJO_38_S6_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SBJO_38_S6_readlength.txt")
SBJO_39_S7_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SBJO_39_S7_readlength.txt")
SBJO_40_S8_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SBJO_40_S8_readlength.txt")
SBJO_41_S9_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SBJO_41_S9_readlength.txt")
SBJO_42_S10_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SBJO_42_S10_readlength.txt")
SJVE_56_S11_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SJVE_56_S11_readlength.txt")
SJVE_57_S12_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SJVE_57_S12_readlength.txt")
SJVE_58_S13_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SJVE_58_S13_readlength.txt")
SJVE_59_S14_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SJVE_59_S14_readlength.txt")
SJVE_60_S15_readlength <- read.table("/storage/group/dut374/default/helaina/data/seq_stats/L_atlanticus_SJVE_60_S15_readlength.txt")


#MODIFY DATAFRAMES
colnames(BBIP_2_S1_readlength) <- "read_length"
colnames(BBIP_3_S2_readlength) <- "read_length"
colnames(BBIP_4_S3_readlength) <- "read_length"
colnames(BBIP_5_S4_readlength) <- "read_length"
colnames(BBIP_6_S5_readlength) <- "read_length"
colnames(SBJO_38_S6_readlength) <- "read_length"
colnames(SBJO_39_S7_readlength) <- "read_length"
colnames(SBJO_40_S8_readlength) <- "read_length"
colnames(SBJO_41_S9_readlength) <- "read_length"
colnames(SBJO_42_S10_readlength) <- "read_length"
colnames(SJVE_56_S11_readlength) <- "read_length"
colnames(SJVE_57_S12_readlength) <- "read_length"
colnames(SJVE_58_S13_readlength) <- "read_length"
colnames(SJVE_59_S14_readlength) <- "read_length"
colnames(SJVE_60_S15_readlength) <- "read_length"


# PLOT DATA
plot_BBIP_2_S1 <- ggplot(BBIP_2_S1_readlength , aes(x = read_length)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_2_S1",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_BBIP_3_S2 <- ggplot(BBIP_3_S2_readlength , aes(x = read_length)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_3_S2",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_BBIP_4_S3 <- ggplot(BBIP_4_S3_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_4_S3",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_BBIP_5_S4 <- ggplot(BBIP_5_S4_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_5_S4",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_BBIP_6_S5 <- ggplot(BBIP_6_S5_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_6_S5",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SBJO_38_S6 <- ggplot(SBJO_38_S6_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_38_S6",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SBJO_39_S7 <- ggplot(SBJO_39_S7_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_39_S7",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SBJO_40_S8 <- ggplot(SBJO_40_S8_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_40_S8",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SBJO_41_S9 <- ggplot(SBJO_41_S9_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_41_S9",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SBJO_42_S10 <- ggplot(SBJO_42_S10_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_42_S10",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SJVE_56_S11 <- ggplot(SJVE_56_S11_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_56_S11",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SJVE_57_S12 <- ggplot(SJVE_57_S12_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_57_S12",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SJVE_58_S13 <- ggplot(SJVE_58_S13_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_58_S13",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SJVE_59_S14 <- ggplot(SJVE_59_S14_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_59_S14",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

plot_SJVE_60_S15 <- ggplot(SJVE_60_S15_readlength, aes(x = V1)) +
  geom_histogram(binwidth = 30, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_60_S15",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(60, 300, by = 30))

# ARRANGE AND SAVE
RL <- ggarrange(
        plot_BBIP_2_S1,
        plot_BBIP_3_S2, 
        plot_BBIP_4_S3, 
        plot_BBIP_5_S4, 
        plot_BBIP_6_S5, 
        plot_SBJO_38_S6, 
        plot_SBJO_39_S7, 
        plot_SBJO_40_S8, 
        plot_SBJO_41_S9, 
        plot_SBJO_42_S10,
        plot_SJVE_56_S11,
        plot_SJVE_57_S12,
        plot_SJVE_58_S13,
        plot_SJVE_59_S14, 
        plot_SJVE_60_S15,
        nrow=5,
        ncol=3)

ggsave("/storage/group/dut374/default/helaina/plots/readlength_distribution.png", RL, width = 12, height = 8)
```

## Run R Script Via Command Line
Do this if your script doesn't need too many resources
```bash
chmod +x plot_readlengths.R
./plot_readlengths.R
```

## Run R Script as Job
```bash
nano plot_readlengths.bash
#!/bin/bash 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=50GB 
#SBATCH --time=3:00:00 
#SBATCH --account=dut374_c
#SBATCH --partition=sla-prio 

#Set Variables
scripts_folder="/storage/group/dut374/default/helaina/scripts"
cd $scripts_folder/plots

#Load R
module use /storage/icds/RISE/sw8/modules/r 
module load r/4.2.1-gcc-8.5.0

#Run R script
R --file=$scripts_folder/plots/plot_readlengths.R
```

## Download
```bash
rsync abc6435@submit.hpc.psu.edu:/storage/home/abc6435/SzpiechLab/abc6435/KROH/plots/readlengths_grid.png /Users/abc6435/Desktop/KROH/plots