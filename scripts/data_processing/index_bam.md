# Index Marked BAM files

## Create Scripts
`samtools index sample.bam sample.bai`
```bash
#Set Variables
scripts_folder="/storage/group/dut374/default/helaina/scripts"

for i in `cat $scripts_folder/ids.txt`; do
    cat<<EOT > $scripts_folder/index_bam_${i}.bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20MB
#SBATCH --time=00:30:00
#SBATCH --account=dut374_c
#SBATCH --job-name=index_bam_${i}
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.err
#SBATCH --output=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

#Set Variables
data_folder="/storage/group/dut374/default/helaina/data"

#loading samtools
module load samtools 

#Index BAM Files
samtools index \$data_folder/bam/${i}_marked.bam \$data_folder/bam/${i}_marked.bai
EOT
done

#Submit each script
ls index_bam* | wc -l 

for i in $scripts_folder/index_bam*.bash; do
    ls ${i}
done

for i in $scripts_folder/index_bam*.bash; do
    sbatch ${i}
done

#Check Job Status
squeue -o "%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R" -u hdl5108

```