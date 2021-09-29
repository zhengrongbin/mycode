#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=03:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=BWAindex # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of cpu cores on one node
#SBATCH --mem=50G

echo "++++ activating biogrids"
source /programs/biogrids.shrc

echo "++++ start runnning for mm10"
bwa index /lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/bwa_index/GRCm39.primary_assembly.genome.fa
echo "++++ start runnning for hg38"
bwa index /lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/hg38/bwa_index/GRCh38.primary_assembly.genome.fa

echo "++++++ finished"

