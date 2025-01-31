# Plot Readlength Distributions

## Load R Module and Install Packages
You only have to install packages once
```bash
module use /storage/icds/RISE/sw8/modules/r 
module load r/4.2.1-gcc-8.5.0

Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('ggplot2', lib='/storage/group/dut374/default/bin/.R')"
```

## Create R Script
```R
scripts_folder="/storage/home/abc6435/SzpiechLab/abc6435/KROH/scripts"
nano $scripts_folder/plots/plot_readlengths.R

#!/usr/bin/env Rscript

# CALL PACKAGES
library(ggplot2, lib.loc=lib='/storage/group/dut374/default/bin/.R')

# IMPORT DATA
ids <- c("L_atlanticus_BBIP_2_S1", "L_atlanticus_BBIP_3_S2", "L_atlanticus_BBIP_4_S3", "L_atlanticus_BBIP_5_S4", "L_atlanticus_BBIP_6_S5",
"L_atlanticus_SBJO_38_S6", "L_atlanticus_SBJO_39_S7", "L_atlanticus_SBJO_40_S8", 
 "L_atlanticus_SBJO_41_S9", "L_atlanticus_SBJO_42_S10",  "L_atlanticus_SJVE_56_S11", "L_atlanticus_SJVE_57_S12", "L_atlanticus_SJVE_58_S13",  "L_atlanticus_SJVE_59_S14", "L_atlanticus_SJVE_60_S15")

base_path <- "/storage/group/dut374/default/helaina/data/seq_stats/"

readlengths <- lapply(ids, function(id) {read.delim(paste0(base_path, id, "_readlength.txt"), header = FALSE)})

###DEBUG THIS ^^ (try doing this in srun)


names(readlengths) <- paste0("readlength_", ids)

# PLOT DATA
RL_183194841 <- ggplot(readlengths$readlength_183194841, aes(x = V1)) +
  geom_histogram(binwidth = 10, fill = "gray") +
  xlab("Read Lengths (bp)") +
  ggtitle(paste("Sample 183194841")) + 
 theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

RL_183195312 <- ggplot(readlengths$readlength_183195312, aes(x = V1)) +
  geom_histogram(binwidth = 10, fill = "gray") +
  xlab("Read Lengths (bp)") +
  ggtitle(paste("Sample 183195312")) + 
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

RL_183195332 <- ggplot(readlengths$readlength_183195332, aes(x = V1)) +
  geom_histogram(binwidth = 10, fill = "gray") +
  xlab("Read Lengths (bp)") +
  ggtitle(paste("Sample 183195332")) + 
  theme_minimal() + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)  
  )

RL_183194861 <- ggplot(readlengths$readlength_183194861, aes(x = V1)) +
  geom_histogram(binwidth = 10, fill = "gray") +
  xlab("Read Lengths (bp)") +
  ggtitle(paste("Sample 183194861")) + 
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

RL_183195321 <- ggplot(readlengths$readlength_183195321, aes(x = V1)) +
  geom_histogram(binwidth = 10, fill = "gray") +
  xlab("Read Lengths (bp)") +
  ggtitle(paste("Sample 183195321")) + 
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

RL_183195304 <- ggplot(readlengths$readlength_183195304, aes(x = V1)) +
  geom_histogram(binwidth = 10, fill = "gray") +
  xlab("Read Lengths (bp)") +
  ggtitle(paste("Sample 183195304")) + 
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

RL_183195326 <- ggplot(readlengths$readlength_183195326, aes(x = V1)) +
  geom_histogram(binwidth = 10, fill = "gray") +
  xlab("Read Lengths (bp)") +
  ggtitle(paste("Sample 183195326")) + 
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black", fill = NA) 
  )

# ARRANGE AND SAVE
RL <- ggarrange(
        RL_183194841,
        RL_183194861, 
        RL_183195304, 
        RL_183195312, 
        RL_183195321, 
        RL_183195326, 
        RL_183195332, 
        nrow=4,
        ncol=3)

ggsave("/storage/home/abc6435/SzpiechLab/abc6435/KROH/plots/readlengths_grid.png", RL, width = 12, height = 8)
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
#SBATCH --account=zps5164_sc 
#SBATCH --partition=sla-prio 

#Set Variables
scripts_folder="/storage/home/abc6435/SzpiechLab/abc6435/KROH/scripts"
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