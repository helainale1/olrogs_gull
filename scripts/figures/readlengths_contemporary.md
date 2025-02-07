# Plot Readlength Distributions

## Load R Module and Install Packages
You only have to install packages once
```bash
module use /storage/icds/RISE/sw8/modules/r 
module load r/4.2.1-gcc-8.5.0

Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('ggplot2', lib='/storage/group/dut374/default/bin/.R')"

Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('scales', lib='/storage/group/dut374/default/bin/.R')"
```

## Create R Script
```R
scripts_folder="/storage/group/dut374/default/helaina/scripts"
nano $scripts_folder/plots/plot_readlengths.R

#!/usr/bin/env Rscript

# CALL PACKAGES
library(ggplot2, lib.loc='/storage/group/dut374/default/bin/.R')

# IMPORT DATA
ids <- c("L_atlanticus_BBIP_2_S1", "L_atlanticus_BBIP_3_S2", "L_atlanticus_BBIP_4_S3", "L_atlanticus_BBIP_5_S4", "L_atlanticus_BBIP_6_S5",
"L_atlanticus_SBJO_38_S6", "L_atlanticus_SBJO_39_S7", "L_atlanticus_SBJO_40_S8", 
 "L_atlanticus_SBJO_41_S9", "L_atlanticus_SBJO_42_S10",  "L_atlanticus_SJVE_56_S11", "L_atlanticus_SJVE_57_S12", "L_atlanticus_SJVE_58_S13",  "L_atlanticus_SJVE_59_S14", "L_atlanticus_SJVE_60_S15")

ids <- c("L_atlanticus_BBIP_2_S1", "L_atlanticus_BBIP_3_S2", "L_atlanticus_BBIP_4_S3")

base_path <- "/storage/group/dut374/default/helaina/data/seq_stats/"

readlengths <- lapply(ids, function(id) {read.delim(paste0(base_path, id, "_readlength.txt"), header = FALSE)})

names(readlengths) <- paste0("readlength_", ids)

# PLOT DATA
RL_L_atlanticus_BBIP_2_S1 <- ggplot(readlengths$readlength_L_atlanticus_BBIP_2_S1, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_2_S1",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_BBIP_3_S2 <- ggplot(readlengths$readlength_L_atlanticus_BBIP_3_S2, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_3_S2",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_BBIP_4_S3 <- ggplot(readlengths$readlength_L_atlanticus_BBIP_4_S3, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_4_S3",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_BBIP_5_S4 <- ggplot(readlengths$readlength_L_atlanticus_BBIP_5_S4, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_5_S4",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_BBIP_6_S5 <- ggplot(readlengths$readlength_L_atlanticus_BBIP_6_S5, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_BBIP_6_S5",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SBJO_38_S6 <- ggplot(readlengths$readlength_L_atlanticus_SBJO_38_S6, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_38_S6",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SBJO_39_S7 <- ggplot(readlengths$readlength_L_atlanticus_SBJO_39_S7, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_39_S7",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SBJO_40_S8 <- ggplot(readlengths$readlength_L_atlanticus_SBJO_40_S8, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_40_S8",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SBJO_41_S9 <- ggplot(readlengths$readlength_L_atlanticus_SBJO_41_S9, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_41_S9",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SBJO_42_S10 <- ggplot(readlengths$readlength_L_atlanticus_SBJO_42_S10, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SBJO_42_S10",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SJVE_56_S11 <- ggplot(readlengths$readlength_L_atlanticus_SJVE_56_S11, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_56_S11",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SJVE_57_S12 <- ggplot(readlengths$readlength_L_atlanticus_SJVE_57_S12, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_57_S12",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SJVE_58_S13 <- ggplot(readlengths$readlength_L_atlanticus_SJVE_58_S13, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_58_S13",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SJVE_59_S14 <- ggplot(readlengths$readlength_L_atlanticus_SJVE_59_S14, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_59_S14",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

RL_L_atlanticus_SJVE_60_S15 <- ggplot(readlengths$readlength_L_atlanticus_SJVE_60_S15, aes(x = V1)) +
  geom_histogram(binwidth = 20, fill = "#8739ea", color = "black") +
  labs(title = "L_atlanticus_SJVE_60_S15",
       x = "Read Length (bp)",
       y = "Frequency") +
scale_x_continuous(breaks = seq(50, 350, by = 50))

# ARRANGE AND SAVE
RL <- ggarrange(
        RL_L_atlanticus_BBIP_2_S1,
        RL_L_atlanticus_BBIP_3_S2, 
        RL_L_atlanticus_BBIP_4_S3, 
        RL_L_atlanticus_BBIP_5_S4, 
        RL_L_atlanticus_BBIP_6_S5, 
        RL_L_atlanticus_SBJO_38_S6, 
        RL_L_atlanticus_SBJO_39_S7, 
        RL_L_atlanticus_SBJO_40_S8, 
        RL_L_atlanticus_SBJO_41_S9, 
        RL_L_atlanticus_SBJO_42_S10,
        RL_L_atlanticus_SJVE_56_S11,
        RL_L_atlanticus_SJVE_57_S12,
        RL_L_atlanticus_SJVE_58_S13,
        RL_L_atlanticus_SJVE_59_S14, 
        RL_L_atlanticus_SJVE_60_S15,
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