# Adapter Removal Pipeline
Using the AdapterRemoval program, we will remove  the adapters at the ends of our reads. The adapters are added to our DNA so that polymerase has a scaffold to start transcription during the PCR step. After transcription, our reads also get the transcribed adapters. But, we don't actually need the adapters and they need to be removed so we can preform our alignments later.  

`path/to/adapter_removal --file1 /path/to/file1 --file2 /path/to/file2 --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /path/to/output/`

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60GB
#SBATCH --time=10:00:00
#SBATCH --account=dut374_c
#SBATCH --partition=sla-prio

set -ue

ID_File="/storage/home/abc6435/ToewsLab/helaina/scripts/ids.txt"
data_folder="/storage/home/abc6435/ToewsLab/helaina/data"
adapter_removal="/storage/home/abc6435/ToewsLab/bin/adapterremoval-2.1.7/build/AdapterRemoval"
	
for  i  in `cat $ID_File`; do echo "$adapter_removal --file1 $data_folder/fastq/${i}_R1_001.fastq.gz --file2 $data_folder/fastq/${i}_R2_001.fastq.gz --collapse --trimns --minlength 20 --qualitybase 33 --gzip --basename $data_folder/trimmed/${i}_trimmed"; done | parallel -j7
```
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20GB
#SBATCH --time=360:00:00
#SBATCH --account=dut374_c
#SBATCH --partition=sla-prio

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_2_S1_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_2_S1_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_BBIP_2_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_3_S2_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_3_S2_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_BBIP_3_trimmed 

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_4_S3_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_4_S3_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_BBIP_4_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_5_S4_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_5_S4_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_BBIP_5_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_6_S5_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_BBIP_6_S5_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_BBIP_6_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_39_S7_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_39_S7_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SBJO_39_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_40_S8_R1_001.fastq.gz --file2 storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_40_S8_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SBJO_40_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_41_S9_R1_001.fastq.gz --file2 storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_41_S9_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SBJO_41_S9_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_42_S10_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SBJO_42_S10_R2_001.fastq.gz --collape --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SBJO_42_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_56_S11_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_56_S11_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SJVE_56_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_57_S12_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_57_S12_R1_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SJVE_57_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_58_S13_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_58_S13_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SJVE_58_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_59_S14_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_59_S14_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SJVE_59_trimmed

/storage/group/dut374/default/bin/adapterremoval-2.1.7/build/AdapterRemoval --file1 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_60_S15_R1_001.fastq.gz --file2 /storage/group/dut374/default/helaina/data/fastq/L_atlanticus_SJVE_60_S15_R2_001.fastq.gz --collapse --trimns --minlength20 --qualitybase 33 --gzip --basename /storage/group/dut374/default/helaina/data/trimmed/L_atlanticus_SJVE_60_trimmed

```
```bash
#Save and write out
^O <ENTER>
^X

#Submit the script to the cluster
sbatch trim_adapters.bash

#To check on your submission
squeue -u hdl5108
```




