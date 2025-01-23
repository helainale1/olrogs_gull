# Converting SAM to BAM

In order to do anything meaningful, we need to convert our SAM files to binary format (BAM). We can accomplish this using SAMTOOLS, which is already pre-installed into the cluster. Specifically, we will use the samtools view command whose function is to convert between human and computer readable formats. We must specify that our input is in SAM format (by default it expects BAM) using the -S option. We must also say that we want the output to be BAM (by default it produces BAM) with the -b option. SAMTOOLS follows the UNIX convention of sending its output to the UNIX STDOUT, so we need to use a redirect operator (“>”) to create a BAM file from the output. http://quinlanlab.org/tutorials/samtools/samtools.html#converting-sam-to-bam-with-samtools-view
samtools view -S -b sample.sam > sample.bam
-S: specifies that the input file is a sam file
-b: species that we want to produce a bam file
>: denotes where to store the output bam file.

```bash 
ID_File="/storage/group/dut374/default/helaina/scripts/ids.txt"
scripts_folder="/storage/group/dut374/default/helaina/scripts"

for i in `cat $ID_File`; do echo "
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --time=5:00:00
#SBATCH --account=dut374_c
#SBATCH --job-name=${i}_create_bam
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.err
#SBATCH --output=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

set -ue
data_folder=\"/storage/group/dut374/default/helaina/data\"

#loading samtools
module load samtools 

#Run samtools
samtools view -S -b \$data_folder/sam/${i}.sam > \$data_folder/bam/${i}.bam" > $scripts_folder/${i}_create_bam.bash; done
```

## Submitting the convert scripts to queue
```bash
#theres an empty line #1 in each script that needs to be removed
for i in `cat $ID_File`; do sed -i '1d' $scripts_folder/${i}_create_bam.bash; done

#submit each script
for i in `cat $ID_File`; do sbatch $scripts_folder/${i}_create_bam.bash; done
```