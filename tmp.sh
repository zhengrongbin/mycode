#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=10:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=m9 # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of cpu cores on one node
#SBATCH --mem=20G


source /programs/biogrids.shrc

echo 'running'

name='MicroC-9_S3'
## merge fastq

cat Fastq/${name}_L001_R1_001.fastq.gz Fastq/${name}_L002_R1_001.fastq.gz Fastq/${name}_L003_R1_001.fastq.gz Fastq/${name}_L004_R1_001.fastq.gz > ${name}.R1.fq.gz
cat Fastq/${name}_L001_R2_001.fastq.gz Fastq/${name}_L002_R2_001.fastq.gz Fastq/${name}_L003_R2_001.fastq.gz Fastq/${name}_L004_R2_001.fastq.gz > ${name}.R2.fq.gz

## run
r1=${name}.R1.fq.gz
r2=${name}.R2.fq.gz

ref='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/bwa_GRCm38/GRCm38.primary_assembly.genome.fa'
genome_file='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/GRCm38.primary_assembly.genome.fa.Genome'
#chromap -i -r $ref -o mm10.chromap.index
## bwa
bwa mem -5SP -T0 -t16 $ref $r1 $r2 > ${name}_bwa.sam

source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate

mkdir temp1
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path $genome_file ${name}_bwa.sam | pairtools sort --tmpdir=./temp1/ --nproc 8 |pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats_${name}.txt|pairtools split --nproc-in 8 --nproc-out 8 --output-pairs mapped_${name}.pairs --output-sam -|samtools view -bS -@16 | samtools sort -@16 -o ${name}_mapped.PT.bam;samtools index ${name}_mapped.PT.bam

jarfile='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20220926_preadip/chromap/juicer_tools_1.22.01.jar'
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/java -Xmx48000m  -Djava.awt.headless=true -jar ${jarfile} pre --threads 8 -q 5 mapped_${name}.pairs microc_${name}.hic mm10.chromSize


echo "Finished"


import pandas as pd
for c in ['TN', 'RT', 'cold2', 'cold7']:
	d = pd.read_csv('/Users/rongbinzheng/Downloads/%s_celltype_avg_exp.csv.gz'%c, compression = 'gzip', index_col = 0)
	d.T.to_csv('%s_celltype_avg_exp.csv.gz'%c, compression = 'gzip')

