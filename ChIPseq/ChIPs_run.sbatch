#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=08:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=CHIPs # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of cpu cores on one node
#SBATCH --mem=50G

echo "++++ activating anaconda"
source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate
conda activate chips

echo "+++++ run CHIPs"

cd /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/ChIPseq

snakemake -s CHIPS/chips.snakefile --rerun-incomplete -j 8 > run.out

echo "+++++ Finished"


