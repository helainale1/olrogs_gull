
# Calculate Coverage with SAMTOOLS
https://www.biostars.org/p/5165/

This will calculate the average depth of coverage across the entire genome. First samtools extracts the depth at every position, sed exludes the mitochondrial genome and any unmapped scaffolds, then awk takes the average depth across all positions. 

## Create Scripts
`samtools depth -a input.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}'`
- `a`: obtains the depth at every position. 
- `awk`: sums the depth at every position and divides by number of positions.

```bash
#Set Variables
scripts_folder="/storage/group/dut374/default/helaina/scripts"
data_folder="/storage/group/dut374/default/helaina/data"

#Make a coverage folder
if [ ! -d "$data_folder/coverage" ]; then
    mkdir -p "$data_folder/coverage"
fi

#Run loop
for i in `cat $scripts_folder/ids.txt`; do
    cat<<EOT > $scripts_folder/calculate_coverage_${i}.bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100MB
#SBATCH --time=05:00:00
#SBATCH --account=dut374_c
#SBATCH --job-name=calculate_coverage_${i}
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.err
#SBATCH --output=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

#Set Variables
data_folder="/storage/group/dut374/default/helaina/data"

#loading samtools
module load samtools 

#Extract Depth
samtools depth -a \$data_folder/bam/${i}_marked.bam > \$data_folder/coverage/${i}_depth.txt

#Remove mito, scafolds, and chrz
sed -i '/^scaffold/d; /^mito/d; /^chrz/d' \$data_folder/coverage/${i}_depth.txt

#Calculate Average Depth
awk '{sum+=\$3} END { print "sample_${i} = " sum/NR }' \$data_folder/coverage/${i}_depth.txt >> \$data_folder/coverage/autosomal_coverage.txt
EOT
done

#Submit Jobs
ls calculate_coverage* | wc -l 

for i in $scripts_folder/calculate_coverage_*.bash; do
     ls ${i}
done

for i in $scripts_folder/calculate_coverage_*.bash; do
    sbatch ${i}
done

#Check Job Status
squeue -o "%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R" -u hdl5108
```