# Call Variants samtools mpileup
Lefouili and Nam 2022 shows the outerperformace of mpileup over GATK HaplotypeCallder after hard filtering. GATK VQSR however was show to outerperform Bcftools mpileup only when coverage is low. Given that I'm calling both contemporary and historical together, I'm going to use mpileup. I'm also using mpileup when calling variants used in GONE, so figure its best to not mix and match tools. 
 

## Alignments
```bash
#Set Variables
scripts_folder="/storage/group/dut374/default/helaina/scripts"
nano $scripts_folder/call_variants_olrogs.bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50GB
#SBATCH --time=80:00:00
#SBATCH --account=dut374_c
#SBATCH --job-name=call_variants_olrogs
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

#load modules
module use /storage/group/dut374/default/sw/modules
module load all

#Set Variables
ref="/storage/group/dut374/default/helaina/data/L_michahellis_ref/bLarMic1.fna"
vcf_dir="/storage/group/dut374/default/helaina/data/vcf"

#Create a bam list
cd /storage/group/dut374/default/helaina/data/bam
realpath *marked.bam > $vcf_dir/bam_list.txt

#Call Variants
bcftools mpileup -f $ref -b $vcf_dir/bam_list.txt | bcftools call -f GQ -mv --ploidy 2 -Oz -o $vcf_dir/olrogs.vcf.gz
```
