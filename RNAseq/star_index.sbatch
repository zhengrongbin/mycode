#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=02:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=STARindex # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of cpu cores on one node
#SBATCH --mem=70G

echo "++++ activating biogrids"
source /programs/biogrids.shrc

echo "++++++ STAR version"
STAR --version

mainpath='/project/RC_Cardio-Chen-e2/rongbinzheng/Genome/mm10'
fapath=${mainpath}/GRCm39.primary_assembly.genome.fa
gtfpath=${mainpath}/gencode.vM27.annotation.gtf

echo "++++ start runnning"

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ${mainpath}/STAR_mm10 \
--genomeFastaFiles ${fapath} --sjdbGTFfile ${gtfpath}
