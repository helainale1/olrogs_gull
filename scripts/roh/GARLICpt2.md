# ROH Calling with GARLIC

## Input Files
```bash
salloc --nodes 1 --ntasks 1 --mem=50G --time=9:00:00 
#Set Variables
data="/storage/group/dut374/default/helaina/data"
vcf="olrogs_tags_bi_qual_nmiss_exhet_renamed_auto.vcf.gz"
garlic="/storage/group/dut374/default/bin/garlic/bin/linux/garlic"
tgls_script="/storage/group/dut374/default/bin/garlic" 

#load modules
module use /storage/group/dut374/default/sw/modules
module load all

#.tped and .tfam
plink --vcf $data/vcf/$vcf --recode transpose --double-id --chr-set 33 --allow-extra-chr --out $data/roh/garlic/olrogs

awk '$1="l_atlanticus"' $data/roh/garlic/olrogs.tfam > $data/roh/garlic/temp && mv -f $data/roh/garlic/temp $data/roh/garlic/olrogs.tfam

#.tgls
$tgls_script/create_tgls_from_vcf.py $data/vcf/$vcf > $data/roh/garlic/olrogs.tgls

#centromere.txt
bcftools query -f '%CHROM\n' $data/vcf/$vcf | uniq > $data/roh/garlic/centromere.txt

awk '$1=$1" 0 1"' $data/roh/garlic/centromere.txt > $data/roh/garlic/temp && mv -f $data/roh/garlic/temp $data/roh/garlic/centromere.txt
```

## Run GARLIC
```bash
#Set Variables
work_dir="/storage/group/dut374/default/helaina/data/roh/garlic"
centromere="/storage/group/dut374/default/helaina/data/roh/garlic/centromere.txt"
garlic_tool="/storage/group/dut374/default/bin/garlic/bin/linux/garlic"

#run 
$garlic_tool --auto-winsize --auto-overlap-frac --winsize 100 --tped $work_dir/olrogs.tped --tfam $work_dir/olrogs.tfam --tgls $work_dir/olrogs.tgls --gl-type GQ --resample 40 --centromere $centromere --size-bounds 1000000 2000000 3000000 4000000 5000000 --out $work_dir/olrogs

#Remove any ROH < 0.5MB
cat $work_dir/olrogs.roh.bed  | awk -F '\t' '$5>500000 || /track/ {print$0}' > $work_dir/olrogs_0.5MB.roh.bed
```

## ROH Summary
```bash
#Set Variables
work_dir="/storage/group/dut374/default/helaina/data/roh/garlic"
shared="/storage/home/abc6435/SzpiechLab/shared"

#Create ID File
sed 's/ /\t/g' $work_dir/olrogs.tfam | cut -f1,2 >> $work_dir/ids.txt
awk '{print $2, $1}' $work_dir/ids.txt > $work_dir/temp && mv -f $work_dir/temp $work_dir/ids.txt

#Calculate ROH Fractions
$tgls_script/calculate_ROH_fractions.pl $work_dir/olrogs.roh.bed 6 $work_dir/ids.txt > $work_dir/olrogs_roh_sum.txt
$tgls_script/calculate_ROH_fractions.pl $work_dir/olrogs_0.5MB.roh.bed 6 $work_dir/ids.txt > $work_dir/olrogs_roh_sum_0.5MB.txt

#Downloading file in local terminal (entered in the local)
rsync hdl5108@submit.hpc.psu.edu:/storage/group/dut374/default/helaina/data/roh/garlic/olrogs_roh_sum_0.5MB.txt /Users/hdl5108/desktop/olrogs_gull/data