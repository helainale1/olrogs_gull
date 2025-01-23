# Marking Duplicated Reads In Sorted BAM Files
Unless you specify, MarkDuplicates does not remove the duplicates, but rather flags them. Removing duplicates would reduce coverage. By flagging them, downstream GATK tools can ignore the reads, and most GATK tools will ignore the duplicates by default.

## Create Scripts
`java -Xmx[ng] -jar picard.jar MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000`

- `-Xmx[ng]`: Memory where n is the number of gigabytes. Must match PBS/SLURM header and if you don't give it sufficient memory, it will kill your job.
- `I`:One or more input SAM or BAM files to analyze. Must be coordinate sorted. Default value: null. This option may be specified 0 or more times.
- `O`: The output file to write marked records to Required.
- `METRICS_FILE`: File to write duplication metrics required.
- `MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000``: Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value: 8000. This option can be set to 'null' to clear the default value.

Note: When I ran this, it took about 25-40 minutes and a total memory of 16GB.

```bash
#Set Variables
scripts_folder="/storage/group/dut374/default/helaina/scripts"
data_folder="/storage/group/dut374/default/helaina/data"

#Run Loop
for i in `cat $scripts_folder/ids.txt`; do
    cat<<EOT > $scripts_folder/mark_duplicates_${i}.bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --time=24:00:00
#SBATCH --account=dut374_c
#SBATCH --job-name=mark_duplicates_${i}
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out
#SBATCH --output=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

#Set Variables
data_folder="/storage/group/dut374/default/helaina/data"
picard_tools_folder="/storage/group/dut374/default/bin/picard_tools_2.20.8"

#Run Picard Tools
java -Xmx20g -jar \$picard_tools_folder/picard.jar MarkDuplicates INPUT=\$data_folder/bam/${i}_sorted.bam OUTPUT=\$data_folder/bam/${i}_marked.bam METRICS_FILE=\$data_folder/bam/${i}_metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000
EOT
done

#Submit each job
for i in $scripts_folder/mark_duplicates_*.bash; do
    ls ${i}
done

for i in $scripts_folder/mark_duplicates_*.bash; do
    sbatch ${i}
done

#Check Job Status
squeue -o "%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R" -u hdl5108
```