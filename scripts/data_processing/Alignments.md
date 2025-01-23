# Aligning the sequences 
## Using bowtie2
Bowtie2 is a good tool for our data because it forms end-to-end read alignments which is best for reads that have already been trimmed for quality and adapters. It is a tool for aligning sequencing reads. 
https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/04_alignment_using_bowtie2.html 
https://bowtie-bio.

## Indexing the Reference Genome
```bash 
nano /storage/group/dut374/default/helaina/scripts/index_ref_genome.bash
#--------------------------------------------NANO---------------------------------
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH --account=dut374_c
#SBATCH --partition=burst  
#SBATCH --qos=burst4x

#Set Variables
gull_ref_folder="/storage/group/dut374/default/helaina/data/L_michahellis_ref"
bowtie2_folder="/storage/group/dut374/default/bin/bowtie2-2.3.5.1"

#Use bowtie to index the reference
$bowtie2_folder/bowtie2-build $gull_ref_folder/GCA_964199755.1_bLarMic1.1_genomic.fna $gull_ref_folder/bLarMic1

#--------------------------------------------EOS---------------------------------
sbatch /storage/group/dut374/default/helaina/scripts/index_ref_genome.bash
squeue -u hdl5108 
```
## Create an alignment script for each sample using a for loop
```bash
for i in `cat /storage/group/dut374/default/helaina/scripts/ids.txt`; do echo "
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200GB
#SBATCH --time=48:00:00
#SBATCH --account=dut374_c
#SBATCH --partition=burst  
#SBATCH --qos=burst4x

#Set Variables
bowtie2_tool=\"/storage/group/dut374/default/bin/bowtie2-2.3.5.1/bowtie2\"
gull_ref=\"/storage/group/dut374/default/helaina/data/L_michahellis_ref/bLarMic1\"
sam_folder=\"/storage/group/dut374/default/helaina/data/sam\"
trimmed_folder=\"/storage/group/dut374/default/helaina/data/trimmed\"
log_folder=\"/storage/group/dut374/default/helaina/job_err_out/bowtie_log\"


#Run Bowtie2
\$bowtie2_tool -p 4 --very-sensitive-local --local -N 0 \
--phred33 -x \$gull_ref --rg-id ${i} --rg SM:${i} \
-1 \$trimmed_folder/${i}_trimmed.pair1.truncated.gz \
-2 \$trimmed_folder/${i}_trimmed.pair2.truncated.gz \
-U \$trimmed_folder/${i}_trimmed.collapsed.gz \
-X 700 -S \$sam_folder/${i}.sam >& \$log_folder/${i}.log" >> /storage/group/dut374/default/helaina/scripts/${i}_align.bash ; done
```

## Use sed and a for loop to remove the black space from each script
```bash
for i in `cat /storage/group/dut374/default/helaina/scripts/ids.txt`; do sed -i '1d' /storage/group/dut374/default/helaina/scripts/${i}_align.bash ; done
```


## Use sbatch and a for loop to submit each script to the cluster
#To Submit the script 
```bash
for i in `cat /storage/group/dut374/default/helaina/scripts/ids.txt`; do sbatch /storage/group/dut374/default/helaina/scripts/${i}_align.bash ; done
```

## To check on your job
```bash
squeue -u hdl5108
```
