## Filtering
I followed a similar filtering protocol as in Wood et al. 2023 for the Bachman's ROH analysis. 

## Testing filters
```bash
nano /storage/group/dut374/default/helaina/scripts/filter_vcf.bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --time=24:00:00
#SBATCH --account=dut374_c
#SBATCH --job-name=filtering_vcf
#SBATCH --error=/storage/group/dut374/default/helaina/job_err_out/%x.%j.err
#SBATCH --output=/storage/group/dut374/default/helaina/job_err_out/%x.%j.out

## Set variables
work_dir="/storage/group/dut374/default/helaina/data/vcf"

#load modules
module use /storage/group/dut374/default/sw/modules
module load all

cd /storage/group/dut374/default/helaina/data/vcf 

#Add All Tags
bcftools +fill-tags -Oz $work_dir/olrogs.vcf.gz -o $work_dir/olrogs_tags.vcf.gz -- -t all

#Biallelic Sites
bcftools view -m2 -M2 -v snps $work_dir/olrogs_tags.vcf.gz -Oz -o $work_dir/olrogs_tags_bi.vcf.gz

#Quality
bcftools filter -e 'QUAL<20' $work_dir/olrogs_tags_bi.vcf.gz -Oz -o $work_dir/olrogs_tags_bi_qual.vcf.gz

#Genotype: Missing
bcftools view -e \
    'N_MISSING>3' \
    $work_dir/olrogs_tags_bi_qual.vcf.gz \
    -Oz -o $work_dir/olrogs_tags_bi_qual_nmiss.vcf.gz

#Excess Heterozygosity (>80%)
bcftools view -i 'COUNT(GT="het")<=12' $work_dir/olrogs_tags_bi_qual_nmiss.vcf.gz -Oz -o $work_dir/olrogs_tags_bi_qual_nmiss_exhet.vcf.gz

#Rename chromosomes
bcftools annotate --rename-chrs $work_dir/chrs.txt -Oz -o $work_dir/olrogs_tags_bi_qual_nmiss_exhet_renamed.vcf.gz $work_dir/olrogs_tags_bi_qual_nmiss_exhet.vcf.gz

#Remove chrZ
bcftools index $work_dir/olrogs_tags_bi_qual_nmiss_exhet_renamed.vcf.gz

awk '{print $2}' $work_dir/chrs.txt | tr '\n' ',' | sed 's/,$//' > chr_list.txt
chrs_list=$(cat chr_list.txt)
bcftools view -r $chrs_list $work_dir/olrogs_tags_bi_qual_nmiss_exhet_renamed.vcf.gz -Oz -o $work_dir/olrogs_tags_bi_qual_nmiss_exhet_renamed_auto.vcf.gz
```

#Genotype: Read Depth and Genotype Quality
bcftools filter \
    -e 'FORMAT/DP < 10 || FORMAT/GQ < 20' \
    --set-GTs . \
    $work_dir/olrogs_tags_bi_qual.vcf.gz \
    -Oz -o $work_dir/olrogs_tags_bi_qual_dp_gq.vcf.gz