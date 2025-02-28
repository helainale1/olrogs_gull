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

awk '{NR=FNR{map[$2]=$1; next}{print map[$1],$0}' chr_names.txt L_atlanticus_BBIP_2_S1_depth.txt
   