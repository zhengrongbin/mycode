#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=03:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=RSEMindex # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of cpu cores on one node
#SBATCH --mem=50G

echo "++++ activating biogrids"
source /programs/biogrids.shrc

echo "++++ running rsem"
rsem-prepare-reference --gtf gencode.vM27.annotation.gtf\
               --bowtie \
               GRCm39.primary_assembly.genome.fa \
               rsem_ref

echo "+++++ Finished"