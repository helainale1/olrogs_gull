# Obtain Read Lengths 
Extract read lengths with SAMTOOLS and then plot the distribution in R. 

## Create Script
```bash
#Set Variables and make a folder
scripts_folder="/storage/group/dut374/default/helaina/scripts"
data_folder="/storage/group/dut374/default/helaina/data"

if [ ! -d "$data_folder/seq_stats" ]; then
    mkdir -p "$data_folder/seq_stats"
fi

nano $scripts_folder/extract_readlenghts.bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem=5GB
#SBATCH --account=dut374_c
#SBATCH --job-name=extract_readlengths
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.err
#SBATCH --output=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

scripts_folder="/storage/group/dut374/default/helaina/scripts"
data_folder="/storage/group/dut374/default/helaina/data"

#loading samtools
module load samtools 

for i in `cat $scripts_folder/ids.txt`; do samtools view $data_folder/bam/${i}.bam | awk '{print length($10)}' >> $data_folder/seq_stats/${i}_readlength.txt; done 
```

## Download to Local Machine
```bash
rsync abc6435@submit.hpc.psu.edu:/storage/home/abc6435/SzpiechLab/abc6435/KROH/data/seq_stats/*.txt /Users/abc6435/Desktop/KROH/data/seq_stats
```

