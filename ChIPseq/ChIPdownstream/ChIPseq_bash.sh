mkdir stats

mv analysis/*/*csv stats/

for i in `ls analysis/align/*/*.sorted.bam|grep -v unique`
do
	rm -rf $i
done


## TWF20
rep1='BAMfiles/TWF20_1_unique.sorted.bam'
rep2='BAMfiles/TWF20_2_unique.sorted.bam'
n='TWF20'

for bam in `ls *bam|grep -v dedup`
do
        echo +++++++$bam
        picard MarkDuplicates I=$bam REMOVE_DUPLICATES=true METRICS_FILE=${bam}_metrices.txt O=${bam}.dedup.bam ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
done


mkdir -p macs/${n}
macs2 callpeak -t ${rep1}.dedup.bam ${rep2}.dedup.bam --SPMR -B -q 0.01 --keep-dup 1 -g hs -n ${n}/${n} -f BAMPE


for i in `ls ../macs/*/*narrowPeak`
do
	j=`basename ${i}`
	awk '{if(($7>5)&&($9>2)){print}}' $i > $j
done
