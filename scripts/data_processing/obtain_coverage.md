```bash
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

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_BBIP_2_S1_depth.txt > temp && mv -f temp L_atlanticus_BBIP_2_S1_depth.txt

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_BBIP_3_S2_depth.txt > temp && mv -f temp L_atlanticus_BBIP_3_S2_depth.txt  

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_BBIP_4_S3_depth.txt > temp && mv -f temp L_atlanticus_BBIP_4_S3_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_BBIP_5_S4_depth.txt > temp && mv -f temp L_atlanticus_BBIP_5_S4_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_BBIP_6_S5_depth.txt > temp && mv -f temp L_atlanticus_BBIP_6_S5_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SBJO_38_S6_depth.txt > temp && mv -f temp L_atlanticus_SBJO_38_S6_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SBJO_39_S7_depth.txt > temp && mv -f temp L_atlanticus_SBJO_39_S7_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SBJO_40_S8_depth.txt > temp && mv -f temp L_atlanticus_SBJO_40_S8_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SBJO_41_S9_depth.txt > temp && mv -f temp L_atlanticus_SBJO_41_S9_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SBJO_42_S10_depth.txt > temp && mv -f temp L_atlanticus_SBJO_42_S10_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SJVE_56_S11_depth.txt > temp && mv -f temp L_atlanticus_SJVE_56_S11_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SJVE_57_S12_depth.txt > temp && mv -f temp L_atlanticus_SJVE_57_S12_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SJVE_58_S13_depth.txt > temp && mv -f temp L_atlanticus_SJVE_58_S13_depth.txt

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SJVE_59_S14_depth.txt > temp && mv -f temp L_atlanticus_SJVE_59_S14_depth.txt 

awk 'NR==FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_SJVE_60_S15_depth.txt > temp && mv -f temp L_atlanticus_SJVE_60_S15_depth.txt 

