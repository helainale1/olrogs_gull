
# Sorting BAM files
When you align FASTQ files with all current sequence aligners, the alignments produced are in random order with respect to their position in the reference genome. In other words, the BAM file is in the order that the sequences occurred in the input FASTQ files. It must be sorted such that the alignments occur in “genome order”. That is, ordered positionally based upon their alignment coordinates on each chromosome.

`samtools sort sample.sam -T sample_temp.bam -o sample_sorted.bam`
- `-T`: Write temporary files. 
- `-o`: Write the final sorted output to the specified file. 

## Create Scripts
```bash
#Set Variables
scripts_folder="/storage/group/dut374/default/helaina/scripts"
data_folder="/storage/group/dut374/default/helaina/data"
ID_File="/storage/group/dut374/default/helaina/scripts/ids.txt"

#Create Loop
for i in `cat $ID_File`; do echo "
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2GB
#SBATCH --time=24:00:00
#SBATCH --account=dut374_c
#SBATCH --job-name=sort_bam_${i}
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.err
#SBATCH --output=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

#Set Variables
data_folder=\"/storage/group/dut374/default/helaina/data\"

#loading samtools
module load samtools

#Run Samtools Sort
samtools sort \$data_folder/bam/${i}.bam -T \$data_folder/bam/${i}_temp.bam -o \$data_folder/bam/${i}_sorted.bam" > $scripts_folder/${i}_sort_bam.bash; done

#theres an empty line #1 in each script that needs to be removed
for i in `cat $ID_File`; do sed -i '1d' $scripts_folder/${i}_sort_bam.bash; done

#Submit each job
#for i in $scripts_folder/sort_bam_*.bash; do
   #sbatch ${i}
#done

#Check Job Status
#squeue -o "%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R" -u hdl5108
```