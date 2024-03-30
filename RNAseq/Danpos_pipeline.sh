bam=$1
Igg=$2
label=$3

source /programs/biogrids.shrc

samtools sort ${bam} --output-fmt BAM > ${bam}.sorted.bam
picard MarkDuplicates -I ${bam}.sorted.bam --REMOVE_DUPLICATES true -M ${bam}_metrices.txt -O ${bam}.dedup.bam ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

dedup_Igg=${Igg}.dedup.bam
if [ -f "$dedup_Igg" ]; then
	echo "+++IgG dedup bam exists!"
else
	samtools sort ${Igg} --output-fmt BAM > ${Igg}.sorted.bam
	picard MarkDuplicates -I ${Igg}.sorted.bam --REMOVE_DUPLICATES true -M ${Igg}_metrices.txt -O ${Igg}.dedup.bam ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
fi

export SAMTOOLS_X=0.1.19

python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${bam}.dedup.bam -b ${dedup_Igg} --smooth_width 0  -o $label -c 25000000 --frsz 200 --extend 200

source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate

chromsize=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize
wigToBigWig -clip ${label}/pooled/${bam}.dedup.bgsub.Fnor.wig $chromsize ${label}/pooled/${bam}.dedup.bgsub.Fnor.wig.bw
