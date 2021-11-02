
### this shell script generate sbatch file that can be submit to process RNA-seq data with paired end sequencing

### metasheet should contains four columns, tab separated
## #1: sample name
## #2 fastq end1 path
## #3 fastq end2 path
## #4 either mm10 or hg38
metasheet=$1  ## a path
while read line
do
    set $line
    name=$1
    path1=$2
    path2=$3
    species=$4

    jobName=$name
    # commpath='/project/RC_Cardio-Chen-e2/'
    commpath='/lab-share/Cardio-Chen-e2/Public/'

    workpath=${commpath}/rongbinzheng/DataProcess/RNA-seq/${jobName}
    QC=${workpath}/QC
    analysis=${workpath}/analysis
    result=${workpath}/result

    version_cont=${analysis}/version.control

    fq1=${path1} #${commpath}/rongbinzheng/BrownAdipo/RNAseq_PreAdipo_Adipo_CL/raw_data/adipocyte+_CL_rep1_1.fq.gz
    fq2=${path2} #${commpath}/rongbinzheng/BrownAdipo/RNAseq_PreAdipo_Adipo_CL/raw_data/adipocyte+_CL_rep1_2.fq.gz

    base_fq1=`basename $fq1`
    base_fq2=`basename $fq2`

    trimed_fq1=${QC}/${base_fq1//.fq.gz/_val_1.fq}
    trimed_fq2=${QC}/${base_fq2//.fq.gz/_val_2.fq}

    transcriptBam=${analysis}/STAR_outputAligned.toTranscriptome.out.bam

    genomicBam=${analysis}/STAR_outputAligned.sortedByCoord.out.bam

    if [ ${species} == hg38 ]
    then
        mainpath=${commpath}/rongbinzheng/Genome/hg38
        fapath=${mainpath}/GRCh38.primary_assembly.genome.fa
        gtfpath=${mainpath}/gencode.v38.annotation.gtf
        indexpath=${mainpath}/STAR_hg38
        salmon_index=${mainpath}/salmon_transcript_index
        transcript_fapath=${mainpath}/gencode.v38.transcripts.fa
    elif [ ${species} == mm10 ]
    then
        mainpath=${commpath}/rongbinzheng/Genome/mm10
        fapath=${mainpath}/GRCm39.primary_assembly.genome.fa
        gtfpath=${mainpath}/gencode.vM27.annotation.gtf
        indexpath=${mainpath}/STAR_mm10
        salmon_index=${mainpath}/salmon_transcript_index
        rsem_ref=${mainpath}/rsem_ref/rsem_ref
        transcript_fapath=${mainpath}/gencode.vM27.transcripts.salmon_quant.fa
    else
        echo '++++++++++ Species Problem +++++++++'
        exit 0
    fi

    echo """#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=04:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=STARrun # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=8 # Number of cpu cores on one node
#SBATCH --mem=70G

echo '++++ activating biogrids'
source /programs/biogrids.shrc

## ==== set input ======
jobName=${jobName}
commpath=${commpath}
fq1=${fq1}
fq2=${fq2}
species=${species}

## ==== set output dir =====
workpath=${commpath}/rongbinzheng/BrownAdipo/analysis/${jobName}
mkdir -p ${workpath}

mkdir -p ${QC}

mkdir -p ${analysis}

echo '++++++ STAR version'
STAR --version

### ========== index path ==========

mainpath=${mainpath}
fapath=${fapath}
gtfpath=${gtfpath}
indexpath=${indexpath}
salmon_index=${salmon_index}
rsem_ref=${rsem_ref}
transcript_fapath=${transcript_fapath}


## =========== process ============
echo '++++ automatic trimming'
trim_galore --version >> ${version_cont}

trim_galore --paired --retain_unpaired --dont_gzip -o $QC --fastqc_args '-d ${QC}' $fq1 $fq2



echo '++++ star runnning'
echo '\n\nSTAR version' >> ${version_cont}
STAR --version >> ${version_cont}

STAR --runThreadN 8 \
--genomeDir ${indexpath} \
--readFilesIn ${trimed_fq1} ${trimed_fq2} \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--quantMode TranscriptomeSAM GeneCounts \
--sjdbGTFfile ${gtfpath} \
--outFileNamePrefix ${analysis}/STAR_output

## ======= counting using RSEM
echo '+++++++ RSEM'
rsem-calculate-expression --version  >> ${version_cont}

rsem-calculate-expression --bam --no-bam-output -p 8 \
    --paired-end --forward-prob 0.5 \
    ${transcriptBam} \
    ${rsem_ref} ${analysis}/rsem

## ======= make bigwig file
echo '++++++++ generating bigwig'
bamCoverage --version >> ${analysis}/version.control
samtools --version|head -1 >> ${analysis}/version.control

samtools index ${genomicBam}

bamCoverage --bam ${genomicBam} --outFileName ${analysis}/${jobName}.bw --normalizeUsing RPKM --outFileFormat bigwig

## ======= change name and move ======
echo '++++++++ move, remove, change result file name'

mkdir -p ${result}
mv ${analysis}/${jobName}.bw ${workpath}/result/${jobName}.bw
mv ${analysis}/rsem.genes.results ${workpath}/result/${jobName}.gene.rsem.txt
mv ${analysis}/STAR_outputReadsPerGene.out.tab ${workpath}/result/${jobName}.gene_count.txt
mv ${analysis}/rsem.isoforms.results ${workpath}/result/${jobName}.isoforms.rsem.txt

rm -rf ${analysis}/STAR_output_STARgenome
rm -rf ${analysis}/STAR_output_STARtmp
rm -rf ${QC}/val*fq
rm -rf ${analysis}/*.toTranscriptome.out.bam 

## ======= annotate genes =====
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/merge_exp_mat.py ann_gene -t ${workpath}/result/${jobName}.gene.rsem.txt -s ${species}
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/merge_exp_mat.py ann_gene -c ${workpath}/result/${jobName}.gene_count.txt -s ${species}

echo '++++++ finished'
""" > ${jobName}.sbatch
done < $metasheet

