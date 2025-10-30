docker pull ubuntu:latest   

curl -H "Accept: application/vnd.github+json" \
 -H "Authorization: zhengrongbin <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>" \
 https://api.github.com/repos/zhengrongbin/MEBOCOST/traffic/clones

chromsize=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize
bam=ChIP_04282023/H2AZ_ab4174_minusCL_unique.sorted.bam
bam=${bam/"/"/_}
label=H2AZ_ab4174_minusCL_vs_IgG

bam2=${bam/"/"/_}
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/wigToBigWig -clip ${label}/pooled/${bam2}.dedup.bgsub.Fnor.wig $chromsize ${label}/pooled/${bam2}.dedup.bgsub.Fnor.wig.bw

## load whole document
docker run -it --name routine -p 8080:8080 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents routine:latest

docker run -it --name seurat -p 8080:8080 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents seurat:v3.1.5

docker run -it --name seurat -p 8081:8081 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents seurat:v2.0.1

## for MECOM
docker run -it --name MECOM -p 8880:8880 -v /Users/rongbinzheng/Documents/BCH/ChenLab:/home/ChenLab -v /Users/rongbinzheng/Documents/CommonData:/home/CommonData seurat:v2.0.1

srun -A cbp -p bch-gpu-pe --qos=normal --gpus=1 --gpus-per-node=1 --mem=30G --pty /bin/bash 
jupyter notebook --no-browser --port 8909 --NotebookApp.iopub_data_rate_limit=10000000000 --ip=0.0.0.0 --NotebookApp.allow_origin=* --allow-root

## for compass
docker run -it --name compass -p 8889:8889 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents  compass:latest
jupyter notebook --allow-root --port=8889 --no-browser --ip=0.0.0.0 

## for vue
docker run -it --name vue -p 8880:8880 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents  vue:latest

jupyter notebook --allow-root --port=8080 --no-browser --ip=0.0.0.0 

remotes::install_version(package = 'Seurat', version = package_version('2.0.1'))

## E2
jupyter notebook --no-browser --port 8900 --NotebookApp.iopub_data_rate_limit=10000000000 --ip=0.0.0.0 --NotebookApp.allow_origin=* --allow-root

ssh -N -L 8900:compute-7-0.rc.tch.harvard.edu:8900 -o ServerAliveInterval=30 ch228298@e3-login.tch.harvard.edu

chr8,HAVANA,gene,83290352,83298452,,+,0,ENSMUSG00000031710.4,protein_coding,Ucp1,2,MGI:98894,OTTMUSG00000061515.1,,,,,,,,,,,

83288352 83292352

awk '{
    if ( $2 == chr8 ) {
        if ( $3 > 83288352 && $3 < 83292352 ) {
            print
        }
    }else if ( $4 == chr8) {
        if ( $5 > 83288352 && $5 < 83292352 ) {
            print
        }
    }
}' chr8.pairs > test.pairs


## mghpcc
## in mghpcc run:
jupyter notebook --no-browser --port 8900 --NotebookApp.iopub_data_rate_limit=10000000000 --ip=0.0.0.0 --NotebookApp.allow_origin=* --allow-root
## in e2 run:
ssh -f -L 8900:compute-m7c1-1-6:8900 ch210487@bch-mghpcc.rc.fas.harvard.edu -N

##
mkdir /lab-share/Cardio-Chen-e2/Public/rongbinzheng/google_drive
rclone mount google-drive:~/ /project/RC_Cardio-Chen-e2/rongbinzheng/google_drive --allow-other --allow-non-empty --vfs-cache-mode writes


srun -A bch -p bch-interactive --pty /bin/bash


### ========== MGHPCC ========== 
## interactive
srun -A bch-mghpcc -p bch-interactive-mghpcc --pty /bin/bash


http://dc2.cistrome.org/batchview/h/44908_45015_44995_62395_3149_3217_3123_3142_3147_44906_47605/w/

### ====== configure homer ========
for i in `ls bin/*.pl`
do
    sed 's/gpfs\/data01\/cbenner\/software/lab-share\/Cardio-Chen-e2\/Public\/rongbinzheng\/anaconda3\/envs\/chips\/share/g' $i > new_bin/$i
done

for i in `ls old/*pm`
do
    sed 's/gpfs\/data01\/cbenner\/software/lab-share\/Cardio-Chen-e2\/Public\/rongbinzheng\/anaconda3\/envs\/chips\/share/g' $i > new_bin/$i
done

snakemake -s CHIPS/chips.snakefile --rerun-incomplete -j 1 -npr


for i in `cat count_file.txt`
do
    python merge_exp_mat.py ann_gene -c $i -s mm10 -d combine/
done

for i in `cat rsem_file.txt`
do
    python merge_exp_mat.py ann_gene -t $i -s mm10 -d combine/
done

python merge_exp_mat.py merge_matrix -i count_file.txt -ft count -d combine/
python merge_exp_mat.py merge_matrix -i rsem_file.txt -ft tpm -d combine/
python merge_exp_mat.py merge_matrix -i rsem_file.txt -ft fpkm -d combine/

## ====
Rscript /Users/rongbinzheng/Documents/CommonCode/RNAseq/diffexp_analysis.R -c count_matrix.csv \
-m comparison_design.csv --gene_ont True --gsea True -s mmu -p ./diffexp/ \
--protein_coding /Users/rongbinzheng/Documents/CommonData/mm10/gencode.vM27.annotation.protein_coding.csv



### gtf sort and zip
perl ../gff3sort/gff3sort.pl --precise --chr_order natural gencode.v19.annotation.gtf | bgzip > hg19.gtf.gz;
tabix -p gff hg19.gtf.gz

### ======= scFEA install
conda create -n scFEA python=3
git clone https://github.com/changwn/scFEA
cd scFEA
conda install --file requirements
conda install pytorch torchvision -c pytorch
pip install --user magic-impute

## test
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file Melissa_full.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv

python src/scFEA.py --data_dir data --input_dir input \
                    --test_file HNSC_GSE103322_expression.csv \
                    --moduleGene_file module_gene_m168.csv \
                    --stoichiometry_matrix cmMat_c70_m168.csv \
                    --sc_imputation True


/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/cold2_mat.csv
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/cold7_mat.csv
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/RT_mat.csv
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/TN_mat.csv


f="TN_mat.csv"
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file $f \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
                    --sc_imputation True 


f="RT_mat.csv"
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file $f \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
                    --sc_imputation True 

f="cold7_mat.csv"
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file $f \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
                    --sc_imputation True 

f="cold2_mat.csv"
python src/scFEA.py --data_dir data --input_dir input \
                    --test_file $f \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
                    --sc_imputation True 



scFEA_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/scFEA/output/RT_mat_module168_cell21153_batch21153_LR0.008_epoch100_SCimpute_T_lambBal1_lambSca1_lambCellCor1_lambModCor_1e-2_20211119-191218.csv
exp_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/RT_mat.csv
meta_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/RT_meta.csv
prefix=BAT_RT

python Metabolic_Comm_Pipeline.py -e $exp_path \
        -m $meta_path \
        -p $scFEA_path \
        -t /lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/test_job/BAT/Metabolic_Communication/data/mouse/metabolite_transporter/manually_check_scFEA.txt \ 
         --cutoff_pvalue 0.05 -d ${prefix} -n ${prefix} --cutoff_exp 0 --cutoff_metabolite 0 --pvalue_method ranksum_pval --cutoff_prop 0.1


scFEA_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/scFEA/output/RT_mat_module168_cell21153_batch21153_LR0.008_epoch100_SCimpute_T_lambBal1_lambSca1_lambCellCor1_lambModCor_1e-2_20211119-191218.csv
exp_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/RT_mat.csv
meta_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/RT_meta.csv
prefix=BAT_TN

scFEA_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/scFEA/output/cold2_mat_module168_cell27078_batch27078_LR0.008_epoch100_SCimpute_T_lambBal1_lambSca1_lambCellCor1_lambModCor_1e-2_20211119-190611.csv
exp_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/cold2_mat.csv
meta_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/cold2_meta.csv
prefix=BAT_cold2

scFEA_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/scFEA/output/cold7_mat_module168_cell19562_batch19562_LR0.008_epoch100_SCimpute_T_lambBal1_lambSca1_lambCellCor1_lambModCor_1e-2_20211119-192337.csv
exp_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/cold7_mat.csv
meta_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/cold7_meta.csv
prefix=BAT_cold7


### cellranger
gsm=GSM4875674
ref_data=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/cellranger/cellranger/mm10/
fastq_path=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/scRNA_BAT/mouse/raw_data/sra

cellranger count --id=${gsm} --transcriptome=${ref_data} --fastqs=${fastq_path}/${gsm} --sample=${gsm} --localcores 1 --r1-length 24



path = 'GSM4875676/GSM4875676_S1_L001_R1_001.fastq.gz'
a_file = gzip.open(path, 'rb')
length = []
for i in a_file:
    length.extend(re.findall('length=[0-9]*', i.decode('utf8').rstrip()))



### build singularity container from docker ones, in e2
module load singularity
singularity build --sandbox compass docker-archive://compass_installed.tar

singularity run /lab-share/Cardio-Chen-e2/Public/rongbinzheng/tmp/compass sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/test_compass/compass_run_meta.sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/test_compass/HNSC_GSE103322_expression.tsv homo_sapiens /lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/test_compass/HNSC_GSE103322_compass_metabolite temp_comp3 14


## merge 
import os,sys
import pandas as pd
path = 'temp_comp3/'
samples = {x:os.path.join(path, x) for x in os.listdir(path) if x.startswith('sample')}
reactions = pd.DataFrame()
for s in samples:
    if not os.path.exists(os.path.join(samples[s], 'uptake.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'uptake.txt'), index_col = 0, sep = '\t')
    reactions = pd.concat([reactions, tmp], axis = 1)

reactions.to_csv('HNSC_compass_uptake_result.csv')

## merge for CCLE compass
import os,sys
import pandas as pd
path = 'ccle_c1_temp'
samples = {x:os.path.join(path, x) for x in os.listdir(path) if x.startswith('sample')}
secretions = pd.DataFrame()
reactions = pd.DataFrame()
uptake = pd.DataFrame()
for s in samples:
    if not os.path.exists(os.path.join(samples[s], 'uptake.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'uptake.txt'), index_col = 0, sep = '\t')
    uptake = pd.concat([uptake, tmp], axis = 1)
    if not os.path.exists(os.path.join(samples[s], 'secretions.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'secretions.txt'), index_col = 0, sep = '\t')
    secretions = pd.concat([secretions, tmp], axis = 1)
    if not os.path.exists(os.path.join(samples[s], 'reactions.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'reactions.txt'), index_col = 0, sep = '\t')
    reactions = pd.concat([reactions, tmp], axis = 1)

path = 'ccle_c2_temp'
samples = {x:os.path.join(path, x) for x in os.listdir(path) if x.startswith('sample')}
for s in samples:
    if not os.path.exists(os.path.join(samples[s], 'uptake.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'uptake.txt'), index_col = 0, sep = '\t')
    uptake = pd.concat([uptake, tmp], axis = 1)
    if not os.path.exists(os.path.join(samples[s], 'secretions.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'secretions.txt'), index_col = 0, sep = '\t')
    secretions = pd.concat([secretions, tmp], axis = 1)
    if not os.path.exists(os.path.join(samples[s], 'reactions.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'reactions.txt'), index_col = 0, sep = '\t')
    reactions = pd.concat([reactions, tmp], axis = 1)


path = 'ccle_c3_temp'
samples = {x:os.path.join(path, x) for x in os.listdir(path) if x.startswith('sample')}
for s in samples:
    if not os.path.exists(os.path.join(samples[s], 'uptake.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'uptake.txt'), index_col = 0, sep = '\t')
    uptake = pd.concat([uptake, tmp], axis = 1)
    if not os.path.exists(os.path.join(samples[s], 'secretions.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'secretions.txt'), index_col = 0, sep = '\t')
    secretions = pd.concat([secretions, tmp], axis = 1)
    if not os.path.exists(os.path.join(samples[s], 'reactions.txt')):
        continue
    tmp = pd.read_csv(os.path.join(samples[s], 'reactions.txt'), index_col = 0, sep = '\t')
    reactions = pd.concat([reactions, tmp], axis = 1)

uptake.to_csv('CCLE_Compass_uptake.csv')
secretions.to_csv('CCLE_Compass_secretions.csv')
reactions.to_csv('CCLE_Compass_reactions.csv')

## to gmt file
out = open('kegg.gmt', 'w')
for line in open('mouse_KEGG_terms_symbol.txt', 'r'):
    line = line.rstrip().split('\t')
    out.write('\t'.join([line[0], 'NA']+line[1].split(';'))+'\n')




((adipose[Title/Abstract]) OR (adipocyte[Title/Abstract]) OR (adipocytes[Title/Abstract]) OR (BAT[Title/Abstract]) OR (WAT[Title/Abstract])) AND ((metabolite sensing[Title/Abstract]) OR (cross-talk[Title/Abstract]) OR (crosstalk[Title/Abstract]) OR (cell communication[Title/Abstract]) OR (cell-cell communication[Title/Abstract]) OR (cellular communication[Title/Abstract]) OR (metabolic communication[Title/Abstract]) OR (metabolite communication[Title/Abstract]) OR (paracrine[Title/Abstract]) OR (autocrine[Title/Abstract]) OR (angiocrine[Title/Abstract]))




samples:
  H2AZ_Ctrl:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/1_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/1_R2_001.fastq.gz
  H2AZ_CL:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/2_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/2_R2_001.fastq.gz
  H2AZac_Ctrl:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/3_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/3_R2_001.fastq.gz
  H2AZac_CL:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/4_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/4_R2_001.fastq.gz
  IgG_Ctrl:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/5_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/5_R2_001.fastq.gz
  IgG_CL:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/6_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/H2AZ_2022/00_fastq/6_R2_001.fastq.gz




fastp 0.20.1
SyntaxError:
Not all output, log and benchmark files of rule peaks_getBroadStats contain the same wildcards. This is crucial though, in order to avoid that two or more jobs write to the same file.
  File "/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/H2AZ_2022_1st/CHIPS/chips.snakefile", line 258, in <module>
  File "/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/H2AZ_2022_1st/CHIPS/modules/peaks.snakefile", line 314, in <module>


./1_R1_001.fastq.gz: OK
md5sum: ./1_R2_001.fastq.gz: No such file or directory
./1_R2_001.fastq.gz: FAILED open or read
md5sum: WARNING: 1 listed file could not be read
md5sum: ./2_R1_001.fastq.gz: No such file or directory
./2_R1_001.fastq.gz: FAILED open or read
md5sum: WARNING: 1 listed file could not be read
md5sum: ./2_R2_001.fastq.gz: No such file or directory
./2_R2_001.fastq.gz: FAILED open or read
md5sum: WARNING: 1 listed file could not be read
md5sum: ./3_R1_001.fastq.gz: No such file or directory
./3_R1_001.fastq.gz: FAILED open or read
md5sum: WARNING: 1 listed file could not be read



rclone sync CommonData google-drive:laptop/CommonData


d = read.table('H2AZ_CL_vs_IgG_pooled_H2AZ_CL.bgsub.Fnor-H2AZ_Ctrl_vs_IgG_pooled_H2AZ_Ctrl.bgsub.Fnor.regions.integrative.xls.heightFC_1.5_P0.05.bed', sep = '\t', header = F)
up = subset(d, V5 > 0)
down = subset(d, V5 < 0)

write.table(up, file = 'H2AZ_CL_vs_Ctrl_IgGSub_dregion_up.bed', sep = '\t', row.names = F, col.names = F, quote = F)
write.table(down, file = 'H2AZ_CL_vs_Ctrl_IgGSub_dregion_down.bed', sep = '\t', row.names = F, col.names = F, quote = F)



python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak T_H:C_H -b T_H:T_IgG,C_H:C_IgG -c T_H:7009569,T_IgG:19124209,C_H:2560138,C_IgG:13591195 --frsz 200 --extend 200 -o TH_vs_CH_dpeak




for i in shcon-BRD4-1.FCHL5FKBBXY_L2_ITTACCGAC-CGAATACG shRB-BRD4-1.FCHL5FKBBXY_L2_ITTCCAGGT-CAGTGCTT shcon-BRD4-2.FCHL5FKBBXY_L2_ITCGTCTGA-GTCCTTGA shRB-BRD4-2.FCHL5FKBBXY_L2_ITACGGTCT-TCCATTGC
do
    bam=${i}.bam
    # mkdir ${i}
    # cd ${i}
    # ln -s /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/CHIPS .
    cp small.sbatch ${i}.sbatch
    echo "samtools fastq -1 ${i}.r1.fq -2 ${i}.r2.fq -0 null -s null -n ../${bam}" >> ${i}.sbatch

done


samples:
   shcon_BRD4_1:
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shcon-BRD4-1.FCHL5FKBBXY_L2_ITTACCGAC-CGAATACG.r1.fq
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shcon-BRD4-1.FCHL5FKBBXY_L2_ITTACCGAC-CGAATACG.r2.fq
  shcon_BRD4_2:
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shcon-BRD4-2.FCHL5FKBBXY_L2_ITCGTCTGA-GTCCTTGA.r1.fq
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shcon-BRD4-2.FCHL5FKBBXY_L2_ITCGTCTGA-GTCCTTGA.r2.fq
  shRB_BRD4_1:
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shRB-BRD4-1.FCHL5FKBBXY_L2_ITTCCAGGT-CAGTGCTT.r1.fq
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shRB-BRD4-1.FCHL5FKBBXY_L2_ITTCCAGGT-CAGTGCTT.r2.fq
  shRB_BRD4_2:
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shRB-BRD4-2.FCHL5FKBBXY_L2_ITACGGTCT-TCCATTGC.r1.fq
   - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/HaojieHuang/shRB-BRD4-2.FCHL5FKBBXY_L2_ITACGGTCT-TCCATTGC.r2.fq


for i in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0
do
    cp noise_small.sbatch noise_${i}.sbatch
    echo "python mebo_cold2_test_noise_replaceTrue.py ${i}" >> noise_${i}.sbatch
    echo "echo 'Finish'" >> noise_${i}.sbatch
done


for i in 0.5 1.5 2 2.5 3
do
    cp stability_small.sbatch stability_${i}.sbatch
    ec=`expr $i*0.866 | bc`
    mc=`expr $i*0.141 | bc`
    echo "python stability_noise_cutoff.py ${ec} ${mc}" >> stability_${i}.sbatch
    echo "echo 'Finish'" >> stability_${i}.sbatch
done


for x in `ls ./|grep -v shELF3|grep T47D`
do
    cat ${x}/sh.sh | sed 's/downgenes/upgenes/g' > ${x}/sh.sh
done

cat KD-minusCL-ATAC_S3_L001_R1_001.fastq.gz KD-minusCL-ATAC_S3_L002_R1_001.fastq.gz \
 KD-minusCL-ATAC_S3_L003_R1_001.fastq.gz KD-minusCL-ATAC_S3_L004_R1_001.fastq.gz > KD_minusCL_ATAC.R1.fq.gz" \
cat KD-minusCL-ATAC_S3_L001_R2_001.fastq.gz KD-minusCL-ATAC_S3_L002_R2_001.fastq.gz \
 KD-minusCL-ATAC_S3_L003_R2_001.fastq.gz KD-minusCL-ATAC_S3_L004_R2_001.fastq.gz > KD_minusCL_ATAC.R2.fq.gz" \

cat KD-plusCL-ATAC_S4_L001_R1_001.fastq.gz KD-plusCL-ATAC_S4_L002_R1_001.fastq.gz \
 KD-plusCL-ATAC_S4_L003_R1_001.fastq.gz KD-plusCL-ATAC_S4_L004_R1_001.fastq.gz > KD_plusCL_ATAC.R1.fq.gz" \
cat KD-plusCL-ATAC_S4_L001_R2_001.fastq.gz KD-plusCL-ATAC_S4_L002_R2_001.fastq.gz \
 KD-plusCL-ATAC_S4_L003_R2_001.fastq.gz KD-plusCL-ATAC_S4_L004_R2_001.fastq.gz > KD_plusCL_ATAC.R2.fq.gz" \

cat shNT-minusCL-ATAC_S1_L001_R1_001.fastq.gz shNT-minusCL-ATAC_S1_L002_R1_001.fastq.gz \
 shNT-minusCL-ATAC_S1_L003_R1_001.fastq.gz shNT-minusCL-ATAC_S1_L004_R1_001.fastq.gz > shNT_minusCL_ATAC.R1.fq.gz" \
cat shNT-minusCL-ATAC_S1_L001_R2_001.fastq.gz shNT-minusCL-ATAC_S1_L002_R2_001.fastq.gz \
 shNT-minusCL-ATAC_S1_L003_R2_001.fastq.gz shNT-minusCL-ATAC_S1_L004_R2_001.fastq.gz > shNT_minusCL_ATAC.R2.fq.gz" \

cat shNT-plusCL-ATAC_S2_L001_R1_001.fastq.gz shNT-plusCL-ATAC_S2_L002_R1_001.fastq.gz \
 shNT-plusCL-ATAC_S2_L003_R1_001.fastq.gz shNT-plusCL-ATAC_S2_L004_R1_001.fastq.gz > shNT_plusCL_ATAC.R1.fq.gz" \
cat shNT-plusCL-ATAC_S2_L001_R2_001.fastq.gz shNT-plusCL-ATAC_S2_L002_R2_001.fastq.gz \
 shNT-plusCL-ATAC_S2_L003_R2_001.fastq.gz shNT-plusCL-ATAC_S2_L004_R2_001.fastq.gz > shNT_plusCL_ATAC.R2.fq.gz" \

samples:
  KD_minusCL_ATAC
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/KD_minusCL_ATAC.R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/KD_minusCL_ATAC.R2.fq.gz" \
  KD_plusCL_ATAC
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/KD_plusCL_ATAC.R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/KD_plusCL_ATAC.R2.fq.gz" \
  shNT_minusCL_ATAC
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/shNT_minusCL_ATAC.R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/shNT_minusCL_ATAC.R2.fq.gz" \
  shNT_plusCL_ATAC
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/shNT_plusCL_ATAC.R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/ATAC/shNT_plusCL_ATAC.R2.fq.gz" \ 


res = li.mt.liana_pipe(adata=adata, groupby = 'bulk_labels',
    resource_name = )

from liana.mt import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

cellchat(adata=adata, groupby = 'bulk_labels')
cellphonedb(adata=adata, groupby = 'bulk_labels')

from liana.mt import rank_aggregate

rank_aggregate(adata = adata, groupby = 'bulk_labels')




h5ad_file = sys.argv[1] ## the path of h5ad file
feature_col = sys.argv[2] ## columns for grouping cell types in adata.obs


adata = sc.read_h5ad(h5ad_file)

if 'process_status' in adata.obs.columns.tolist():
    adata = adata[adata.obs['process_status'] == 'QC pass']

adata.obs['disease'] = ['Normal' if x == 'NA' else x for x in adata.obs['disease'].tolist()]

diseases = adata.obs['disease'].unique().tolist()
for d in diseases:
    file = 'liana_result/%s.%s.liana.csv'%(os.path.basename(h5ad_file), d.replace('_', '').replace(' ', '_').replace(',', '').replace("'", '').replace("’", '').replace('/', '_or_'))
    if os.path.exists(file):
        continue
    ### 1. pass expression data by scanpy adata object
    adata_new = adata[adata.obs['disease'] == d]
    ct_count = adata_new.obs.groupby(['disease', 'ct'])['orig.ident'].count()
    if ct_count[ct_count>5].shape[0] == 1:
        continue
    adata_new.obs = adata_new.obs[[feature_col]]
    adata_new.obs[feature_col] = pd.Categorical(adata_new.obs['ct'])
    res = rank_aggregate(adata = adata_new, groupby = feature_col,
                         use_raw = False, inplace = False)
    res.to_csv(file)


for i in adrenal_gland bone_marrow brain breast breast_milk eye gut heart kidney liver lung pancreas skin thymus
do
    cp simple.sbatch ${i}_liana.sbatch
    echo 'python liana_run.py data/'${i}'.h5ad ct' >> ${i}_liana.sbatch
    echo "echo Finished" >> ${i}_liana.sbatch
done




d = sc.read_h5ad('data/adipose.h5ad')
d1 = d.obs[d.obs['project_id'] == 'GSE129363'][['sample', 'sample_type', 'tissue', 'disease', 'anatomical_site', 'bmi', 'ct']]
d1[~d1['sample'].duplicated()]
d1[(d1['anatomical_site'] == 'visceral adipose tissue') & (d1['sample_type'] == 'normal')]

d2 = d.obs[d.obs['project_id'] == 'GSE136229'][['sample', 'sample_type', 'tissue', 'disease', 'anatomical_site', 'bmi', 'ct']]
d2[~d2['sample'].duplicated()]

dd = d.obs[~d.obs['sample'].duplicated()][['sample', 'project_id', 'sample_type', 'tissue', 'disease', 'anatomical_site', 'bmi', 'ct', 'cell_sorting']].sort_values('sample')

                                    sample project_id               sample_type   tissue   disease              anatomical_site    bmi                        ct          cell_sorting
AAACCTGAGCCATCGC-1--GSM3711757  GSM3711757  GSE129363  non-tumor disease tissue  adipose  diabetes  subcutaneous adipose tissue  obese                  pericyte                    NA
AAACCTGAGATGCCAG-1--GSM3711758  GSM3711758  GSE129363  non-tumor disease tissue  adipose  diabetes      visceral adipose tissue  obese             CD14 monocyte                    NA
AAACCTGAGCCACTAT-1--GSM3711759  GSM3711759  GSE129363  non-tumor disease tissue  adipose  diabetes  subcutaneous adipose tissue  obese                macrophage                    NA
AAACCTGAGAAAGTGG-1--GSM3711760  GSM3711760  GSE129363  non-tumor disease tissue  adipose  diabetes      visceral adipose tissue  obese  adipose ITLN1 fibroblast                    NA
AAACCTGCACCGCTAG-1--GSM3711771  GSM3711771  GSE129363  non-tumor disease tissue  adipose  diabetes      visceral adipose tissue  obese  adipose ITLN1 fibroblast                    NA
AAACCTGCAATCAGAA-1--GSM3711773  GSM3711773  GSE129363                    normal  adipose        NA      visceral adipose tissue  obese  adipose ITLN1 fibroblast                    NA
AAACCTGTCTGTCCGT-1--GSM3711774  GSM3711774  GSE129363                    normal  adipose        NA  subcutaneous adipose tissue  obese                   CD16 NK                    NA
AAACCTGCACCTCGTT-1--GSM3711775  GSM3711775  GSE129363                    normal  adipose        NA      visceral adipose tissue  obese                   CD16 NK                    NA
AAACGGGAGTGTTAGA-1--GSM3711776  GSM3711776  GSE129363  non-tumor disease tissue  adipose  diabetes      visceral adipose tissue  obese                GZMK CD8 T                    NA
AAACCTGGTAAACCTC-1--GSM3711777  GSM3711777  GSE129363  non-tumor disease tissue  adipose  diabetes  subcutaneous adipose tissue  obese                 adipocyte                    NA
AAACGGGAGCTGTCTA-1--GSM3711779  GSM3711779  GSE129363                    normal  adipose        NA      visceral adipose tissue  obese  adipose ITLN1 fibroblast                    NA
AAACCTGCACAGACTT-1--GSM3711780  GSM3711780  GSE129363                    normal  adipose        NA      visceral adipose tissue  obese  adipose ITLN1 fibroblast                    NA
AAACCTGCAGTCTTCC-1--GSM3711781  GSM3711781  GSE129363                    normal  adipose        NA  subcutaneous adipose tissue  obese                    mregDC                    NA
AAACCTGCATAACCTG-1--GSM3711782  GSM3711782  GSE129363                    normal  adipose        NA      visceral adipose tissue  obese                GZMK CD8 T                 CD34-
AAACCTGAGATGTGTA-1--GSM3711783  GSM3711783  GSE129363                    normal  adipose        NA      visceral adipose tissue  obese            CFD fibroblast                 CD34+
AAACCTGAGCACACAG-1--GSM3711784  GSM3711784  GSE129363                    normal  adipose        NA  subcutaneous adipose tissue  obese            CFD fibroblast                 CD34-
AAACCTGAGCTCAACT-1--GSM3711786  GSM3711786  GSE129363  non-tumor disease tissue  adipose  diabetes      visceral adipose tissue  obese               naive CD4 T                 CD34-
AAACCTGAGAGTACCG-1--GSM3711787  GSM3711787  GSE129363  non-tumor disease tissue  adipose  diabetes      visceral adipose tissue  obese            CFD fibroblast                 CD34+
AAACCTGAGACACGAC-1--GSM3711788  GSM3711788  GSE129363  non-tumor disease tissue  adipose  diabetes  subcutaneous adipose tissue  obese            CFD fibroblast                 CD34+
AAACCTGAGTTGTAGA-1--GSM4042945  GSM4042945  GSE136229                    normal  adipose        NA      visceral adipose tissue     NA  adipose ITLN1 fibroblast                    NA
AAACGGGAGCAGCGTA-1--GSM4042946  GSM4042946  GSE136229                    normal  adipose        NA      visceral adipose tissue     NA  adipose ITLN1 fibroblast                    NA
AAACCCAAGCGAAACC-1--GSM4717158  GSM4717158  GSE155960                    normal  adipose        NA                           NA     NA     KLRB1 cytotoxic CD4 T                 CD45+
AAACCCAAGGTTTACC-1--GSM4717159  GSM4717159  GSE155960                    normal  adipose        NA                           NA     NA                GZMK CD8 T                 CD45+
AAACCCATCAAGTTGC-1--GSM4717160  GSM4717160  GSE155960                    normal  adipose        NA                           NA     NA              memory CD4 T                 CD45+
AAACCCAAGCATTGTC-1--GSM4717161  GSM4717161  GSE155960                    normal  adipose        NA                           NA     NA              memory CD4 T                 CD45+
AAACCCAGTGTTTACG-1--GSM4717162  GSM4717162  GSE155960                    normal  adipose        NA                           NA     NA              memory CD4 T                 CD45+
AAACCCAAGACAGCGT-1--GSM4717163  GSM4717163  GSE155960                    normal  adipose        NA                           NA     NA              memory CD4 T                 CD45+
AAACCCAAGTTACGAA-1--GSM4724861  GSM4724861  GSE156110                    normal  adipose        NA                           NA     NA                   CD16 NK   Natural killer cell
AAACCCACAGAGGAAA-1--GSM4724862  GSM4724862  GSE156110                    normal  adipose        NA                           NA     NA                   CD16 NK   Natural killer cell
AAACCCAAGCACTTTG-1--GSM4724863  GSM4724863  GSE156110                    normal  adipose        NA                           NA     NA                   CD16 NK  Innate lymphoid cell
AAACCCACAATCGAAA-1--GSM4724864  GSM4724864  GSE156110                    normal  adipose        NA                           NA     NA           SLAMF1-XCL+ ILC  Innate lymphoid cell
AAACCCAAGAGCAAGA-1--GSM4724865  GSM4724865  GSE156110                    normal  adipose        NA                           NA     NA                      cDC2        Dendritic cell
AAACCCACACCCAAGC-1--GSM4724866  GSM4724866  GSE156110                    normal  adipose        NA                           NA     NA                      cDC2        Dendritic cell
AAACGAATCTGTCCGT-1--GSM4724867  GSM4724867  GSE156110                    normal  adipose        NA                           NA     NA             CD16 monocyte            Macrophage
AAACCCATCTAGGCCG-1--GSM4724868  GSM4724868  GSE156110                    normal  adipose        NA                           NA     NA             CD14 monocyte            Macrophage
AAACCCATCAAGTCGT-1--GSM5436518  GSM5436518  GSE179887               solid tumor  adipose    lipoma    stromal vascular fraction     NA             CD14 monocyte                    NA
AAACCCAAGAAATTGC-1--GSM5436519  GSM5436519  GSE179887               solid tumor  adipose    lipoma    stromal vascular fraction     NA                GZMB CD8 T                    NA
>>> 



/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20221226_MicroC_Yang/raw_data/12.17.22.MicroC/shNT-minusCL-rep1_S1_mapped.PT.bam
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20221226_MicroC_Yang/raw_data/12.17.22.MicroC/shNT-minusCL-rep2_S2_mapped.PT.bam
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20230106_MicroC_Yang/MicroC-9_S1_mapped.PT.bam
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20230113_MicroC_Yang/MicroC-13_S1_mapped.PT.bam


/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20221226_MicroC_Yang/raw_data/12.26.22.MicroC/MicroC-2_S1_mapped.PT.bam
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20221226_MicroC_Yang/raw_data/12.26.22.MicroC/MicroC-6_S2_mapped.PT.bam
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20230107_MicroC_Yang/MicroC-10_S3_mapped.PT.bam
/lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/MicroC/20230113_MicroC_Yang/MicroC-14_S2_mapped.PT.bam



run-merge-pairs.sh shNT_minusCL_merge_mapq40_pairs ./pairs/mapped_shNT-minusCL-rep1_S1.pairs.gz ./pairs/mapped_shNT-minusCL-rep2_S2.pairs.gz ./pairs/mapped_MicroC-9_S1.pairs.gz ./pairs/mapped_MicroC-13_S1.pairs.gz 

run-merge-pairs.sh shNT_plusCL_merge_mapq40_pairs ./pairs/mapped_MicroC-2_S1.pairs.gz ./pairs/mapped_MicroC-6_S2.pairs.gz ./pairs/mapped_MicroC-10_S3.pairs.gz ./pairs/mapped_MicroC-14_S2.pairs.gz


STRIDE deconvolve --sc-count ../sc_counts.tsv --sc-celltype ../sc_meta.tsv --st-count ../st_counts.tsv --normalize

STRIDE deconvolve --sc-count ../sc_count.tsv --sc-celltype ../sc_meta.tsv --st-count ../st_count.tsv --normalize

STRIDE deconvolve --sc-count ../sc_count_hum.tsv --sc-celltype ../sc_meta.tsv --st-count ../st_count_hum.tsv --normalize

cat Fastq/RNA${i}_S${i}_L001_R1_001.fastq.gz Fastq/RNA${i}_S${i}_L002_R1_001.fastq.gz Fastq/RNA${i}_S${i}_L003_R1_001.fastq.gz Fastq/RNA${i}_S${i}_L004_R1_001.fastq.gz > RNA${i}_R1.fq.gz" \
cat Fastq/RNA${i}_S${i}_L001_R2_001.fastq.gz Fastq/RNA${i}_S${i}_L002_R2_001.fastq.gz Fastq/RNA${i}_S${i}_L003_R2_001.fastq.gz Fastq/RNA${i}_S${i}_L004_R2_001.fastq.gz > RNA${i}_R2.fq.gz" \


samples:
  H2AZ_minusCL_rep1:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-1_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-1_R2.fq.gz" \
  IgG_minusCL_rep1:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-2_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-2_R2.fq.gz" \
  H2AZ_plusCL_rep1:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-3_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-3_R2.fq.gz" \
  IgG_plusCL_rep1:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-4_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-4_R2.fq.gz" \
  H2AZ_minusCL_rep2:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-5_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-5_R2.fq.gz" \
  IgG_minusCL_rep2:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-6_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-6_R2.fq.gz" \
  H2AZ_plusCL_rep2:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-7_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-7_R2.fq.gz" \
  IgG_plusCL_rep2:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-8_R1.fq.gz" \
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/H2AZ_MNase/DataProcess/CUTRUN/Fastq/CUTandRUN-8_R2.fq.gz" \

for i in MicroC-13_S1 MicroC-9_S1 shNT-minusCL-rep1_S1 shNT-minusCL-rep2_S2
do
    fq1=Fastq/minusCL/${i}.R1.fq.gz" \
    fq2=Fastq/minusCL/${i}.R2.fq.gz" \
    mkdir -p $i
    trim_galore --paired -o $i $fq1 $fq2
done

for i in MicroC-10_S3 MicroC-14_S2 MicroC-2_S1 MicroC-6_S2
do
    fq1=Fastq/plusCL/${i}.R1.fq.gz" \
    fq2=Fastq/plusCL/${i}.R2.fq.gz" \
    mkdir -p $i
    trim_galore --paired -o $i $fq1 $fq2
done


mageck count -l ../../H3_design_addgene.txt --fastq "/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/Nalm6_sample_2/Nalm6_sample_2_USPD16093943-AK4836_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/Nalm6_sample_3/Nalm6_sample_3_USPD16093943-AK4837_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/Nalm6_sample_4/Nalm6_sample_4_USPD16093943-AK74_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/Nalm6_sample_5/Nalm6_sample_5_USPD16093943-AK4838_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/Nalm6_sample_6/Nalm6_sample_6_USPD16093943-AK4839_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/Nalm6_sample_13/Nalm6_sample_13_USPD16093943-AK4951_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/Nalm6_sample_14/Nalm6_sample_14_USPD16093943-AK4952_HWNJJBBXX_L7_1.fq.gz" \
--sample-label sample2,sample3,sample4,sample5,sample6,sample13,sample14 --norm-method control --control-sgrna ../../AAVS1.txt -n count/Nalm6_121218


mageck count -l ../../H3_design_addgene.txt --fastq "/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/s697_sample_7/s697_sample_7_USPD16093943-AK4840_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/s697_sample_8/s697_sample_8_USPD16093943-AK4841_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/s697_sample_9/s697_sample_9_USPD16093943-AK4842_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/s697_sample_10/s697_sample_10_USPD16093943-AK4843_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/s697_sample_11/s697_sample_11_USPD16093943-AK4844_HWNJJBBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_121218/s697_sample_12/s697_sample_12_USPD16093943-AK4845_HWNJJBBXX_L7_1.fq.gz" \
--sample-label sample7,sample8,sample9,sample10,sample11,sample12 --norm-method control --control-sgrna ../../AAVS1.txt -n count/s697_121218


mageck count -l ../../H3_design_addgene.txt --fastq "/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/Nalm6_1/Nalm6_1_USPD16094204-AK7408_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/Nalm6_2/Nalm6_2_USPD16094204-AK4836_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/Nalm6_3/Nalm6_3_USPD16094204-AK4837_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/Nalm6_4/Nalm6_4_USPD16094204-AK74_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/Nalm6_5/Nalm6_5_USPD16094204-AK4838_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/Nalm6_6/Nalm6_6_USPD16094204-AK4839_HWKH5BBXX_L7_1.fq.gz" \
--sample-label sample1,sample2,sample3,sample4,sample5,sample6 --norm-method control --control-sgrna ../../AAVS1.txt -n count/Nalm6_122418

mageck count -l ../../H3_design_addgene.txt --fastq "/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/s697_7/s697_7_USPD16094204-AK4840_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/s697_8/s697_8_USPD16094204-AK4841_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/s697_9/s697_9_USPD16094204-AK4842_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/s697_10/s697_10_USPD16094204-AK4843_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/s697_11/s697_11_USPD16094204-AK4844_HWKH5BBXX_L7_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Qin/screen_122418/s697_12/s697_12_USPD16094204-AK4845_HWKH5BBXX_L7_1.fq.gz" \
--sample-label sample7,sample8,sample9,sample10,sample11,sample12 --norm-method control --control-sgrna ../../AAVS1.txt -n count/s697_122418


mageck count -l ../library.txt --fastq "/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/RCH-ACV/DMSO1_USPD16092791_HW2M2CCXY_L1_1.fq.gz,/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/RCH-ACV/DMSO2_USPD16092797_HW2M2CCXY_L2_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/RCH-ACV/P9B_5D1_USPD16092789_HW2M2CCXY_L1_1.fq.gz,/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/RCH-ACV/P9B_5D2_USPD16092795_HW2M2CCXY_L2_1.fq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/RCH-ACV/dex1_USPD16092792_HW2M2CCXY_L1_1.fq.gz,/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/RCH-ACV/dex2_USPD16092798_HW2M2CCXY_L2_1.fq.gz" \
--sample-label DMSO,P9B_5D,Dex --norm-method control --control-sgrna ../non-target.txt -n count/RCHACV

mageck count -l ../library.txt --fastq "/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/SEM/SEM-DMSO-1_S6_L001_R1_001.fastq.gz,/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/SEM/SEM-DMSO-2_S7_L001_R1_001.fastq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/SEM/SEM-Dex-1_S8_L001_R1_001.fastq.gz,/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/SEM/SEM-Dex-2_S9_L001_R1_001.fastq.gz" \
"/Volumes/GoogleDrive/My Drive/Collaboration/BcellALL/Licht/SEM/SEM-Ref_S5_L001_R1_001.fastq.gz" \
--sample-label DMSO,Dex,Ref --norm-method control --control-sgrna ../non-target.txt -n count/SEM

for i in 5_S2 6_S3 7_S4 8_S5
do
cat 01.07.2023.Fastq/ATAC-${i}_L001_R1_001.fastq.gz 01.07.2023.Fastq/ATAC-${i}_L002_R1_001.fastq.gz 01.07.2023.Fastq/ATAC-${i}_L003_R1_001.fastq.gz 01.07.2023.Fastq/ATAC-${i}_L004_R1_001.fastq.gz > ATAC${i}_R1.fq.gz
cat 01.07.2023.Fastq/ATAC-${i}_L001_R2_001.fastq.gz 01.07.2023.Fastq/ATAC-${i}_L002_R2_001.fastq.gz 01.07.2023.Fastq/ATAC-${i}_L003_R2_001.fastq.gz 01.07.2023.Fastq/ATAC-${i}_L004_R2_001.fastq.gz > ATAC${i}_R2.fq.gz
done

shNT_minusCL_rep1
  - /temp_work/ch228298/ATAC/ATAC1_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC1_R2.fq.gz
shNT_plusCL_rep1
  - /temp_work/ch228298/ATAC/ATAC2_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC2_R2.fq.gz
KD_minusCL_rep1
  - /temp_work/ch228298/ATAC/ATAC3_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC3_R2.fq.gz
KD_plusCL_rep1
  - /temp_work/ch228298/ATAC/ATAC4_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC4_R2.fq.gz
shNT_minusCL_rep2
  - /temp_work/ch228298/ATAC/ATAC5_S2_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC5_S2_R2.fq.gz
shNT_plusCL_rep2
  - /temp_work/ch228298/ATAC/ATAC6_S3_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC6_S3_R2.fq.gz
KD_minusCL_rep2
  - /temp_work/ch228298/ATAC/ATAC7_S4_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC7_S4_R2.fq.gz
KD_plusCL_rep2
  - /temp_work/ch228298/ATAC/ATAC8_S5_R1.fq.gz
  - /temp_work/ch228298/ATAC/ATAC8_S5_R2.fq.gz

RunName,Treat1,Cont1,Treat2,Cont2
shNT_minusCL_rep1,shNT_minusCL_rep1,,,
shNT_plusCL_rep1,shNT_plusCL_rep1,,,
KD_minusCL_rep1,KD_minusCL_rep1,,,
KD_plusCL_rep1,KD_plusCL_rep1,,,
shNT_minusCL_rep2,shNT_minusCL_rep2,,,
shNT_plusCL_rep2,shNT_plusCL_rep2,,,
KD_minusCL_rep2,KD_minusCL_rep2,,,
KD_plusCL_rep2,KD_plusCL_rep2,,,


Samples baseline    DMSO    Dex
Ref
DMSO
Dex

mageck mle -k ./count/SEM.count.txt -d design.txt --control-sgrna ../non-target.txt --norm-method control --adjust-method fdr -n SEM_mle


flux_res = {} #pd.DataFrame()
for s in os.listdir(tmp):
    ss = pd.read_csv(tmp+'/'+s, index_col = 0)
    ss.columns = [s.replace('.csv', '')]
    # flux_res = pd.concat([flux_res, ss], axis = 1)
    flux_res[s] = ss


docker run --detach \
           --publish 8989:80 \
           --volume /Users/rongbinzheng/Documents/TsengLab/H2AZ_Yang/MicroC-2023/mcool:/data \
           --name higlass-container \
         higlass/higlass-docker:v0.6.1

docker exec higlass-container python higlass-server/manage.py \
  ingest_tileset \
  --filename /data/microc_shNT_minusCL_merge.from_hic.mcool \
  --datatype matrix \
  --filetype cooler 


/lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/java -Xmx10g -jar ${jarfile} hiccupsdiff --cpu -k KR microc_shNT_plusCL_merge_mapq5.hic microc_shNT_minusCL_merge_mapq5.hic microc_shNT_plusCL_merge_mapq5_hiccups/merged_loops.bedpe microc_shNT_minusCL_merge_mapq5_hiccups/merged_loops.bedpe diff_hiccups_plusCL_vs_minusCL

python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/mustache/mustache/diff_mustache.py -f1 microc_shNT_plusCL_merge_mapq5.hic -f2 microc_shNT_minusCL_merge_mapq5.hic -r 5000 -st 0.8 -o diff_mustache_plusCL_vs_minusCL -pt 0.05 -pt 0.1 -p 8

fanc compartments -f -v microc_shNT_minusCL_merge_mapq5.from_hic.mcool.10kb_egvector_fanc.bed -d microc_shNT_minusCL_merge_mapq5.from_hic.mcool.10kb_comp_fanc.bed -g GRCm38.primary_assembly.genome.fa microc_shNT_minusCL_merge_mapq5.from_hic.mcool@10000 microc_shNT_minusCL_merge_mapq5.from_hic.mcool.10kb_comp_fanc.ab

cut -f 1-3 diff_mustache_plusCL_vs_minusCL.loop1|grep -v 'BIN' > loop1.anchor.bed  
cut -f 4-6 diff_mustache_plusCL_vs_minusCL.loop1|grep -v 'BIN'  >> loop1.anchor.bed

cut -f 1-3 diff_mustache_plusCL_vs_minusCL.loop2|grep -v 'BIN' > loop2.anchor.bed  
cut -f 4-6 diff_mustache_plusCL_vs_minusCL.loop2|grep -v 'BIN'  >> loop2.anchor.bed

cut -f 1-3 diff_mustache_plusCL_vs_minusCL.diffloop2|grep -v 'BIN' > diffloop2.anchor.bed  
cut -f 4-6 diff_mustache_plusCL_vs_minusCL.diffloop2|grep -v 'BIN'  >> diffloop2.anchor.bed

cat loop1.anchor.bed loop2.anchor.bed diffloop1.anchor.bed diffloop2.anchor.bed|sort|uniq > all_loops.bed


cat Fastq/CUTRUN-1_S3_L001_R1_001.fastq.gz Fastq/CUTRUN-1_S3_L002_R1_001.fastq.gz Fastq/CUTRUN-1_S3_L003_R1_001.fastq.gz Fastq/CUTRUN-1_S3_L004_R1_001.fastq.gz > Fastq/H2AZ_ActiveMotif_0C.R1.fq.gz
cat Fastq/CUTRUN-1_S3_L001_R2_001.fastq.gz Fastq/CUTRUN-1_S3_L002_R2_001.fastq.gz Fastq/CUTRUN-1_S3_L003_R2_001.fastq.gz Fastq/CUTRUN-1_S3_L004_R2_001.fastq.gz > Fastq/H2AZ_ActiveMotif_0C.R2.fq.gz

cat Fastq/CUTRUN-2_S4_L001_R1_001.fastq.gz Fastq/CUTRUN-2_S4_L002_R1_001.fastq.gz Fastq/CUTRUN-2_S4_L003_R1_001.fastq.gz Fastq/CUTRUN-2_S4_L004_R1_001.fastq.gz > Fastq/IgG_ActiveMotif_0C.R1.fq.gz
cat Fastq/CUTRUN-2_S4_L001_R2_001.fastq.gz Fastq/CUTRUN-2_S4_L002_R2_001.fastq.gz Fastq/CUTRUN-2_S4_L003_R2_001.fastq.gz Fastq/CUTRUN-2_S4_L004_R2_001.fastq.gz > Fastq/IgG_ActiveMotif_0C.R2.fq.gz

cat Fastq/CUTRUN-3_S5_L001_R1_001.fastq.gz Fastq/CUTRUN-3_S5_L002_R1_001.fastq.gz Fastq/CUTRUN-3_S5_L003_R1_001.fastq.gz Fastq/CUTRUN-3_S5_L004_R1_001.fastq.gz > Fastq/H2AZ_ActiveMotif_4C.R1.fq.gz
cat Fastq/CUTRUN-3_S5_L001_R2_001.fastq.gz Fastq/CUTRUN-3_S5_L002_R2_001.fastq.gz Fastq/CUTRUN-3_S5_L003_R2_001.fastq.gz Fastq/CUTRUN-3_S5_L004_R2_001.fastq.gz > Fastq/H2AZ_ActiveMotif_4C.R2.fq.gz

cat Fastq/CUTRUN-4_S6_L001_R1_001.fastq.gz Fastq/CUTRUN-4_S6_L002_R1_001.fastq.gz Fastq/CUTRUN-4_S6_L003_R1_001.fastq.gz Fastq/CUTRUN-4_S6_L004_R1_001.fastq.gz > Fastq/IgG_ActiveMotif_4C.R1.fq.gz
cat Fastq/CUTRUN-4_S6_L001_R2_001.fastq.gz Fastq/CUTRUN-4_S6_L002_R2_001.fastq.gz Fastq/CUTRUN-4_S6_L003_R2_001.fastq.gz Fastq/CUTRUN-4_S6_L004_R2_001.fastq.gz > Fastq/IgG_ActiveMotif_4C.R2.fq.gz


  H2AZ_ab4174:
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN1_S3.R1.fq.gz
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN1_S3.R2.fq.gz
  H2AZ_ab188314:
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN2_S4.R1.fq.gz
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN2_S4.R2.fq.gz
  H2AZ_ab124793:
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN3_S5.R1.fq.gz
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN3_S5.R2.fq.gz
  H2AZ_ActiveMotif:
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN4_S6.R1.fq.gz
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN4_S6.R2.fq.gz
  H3K4me3:
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN5_S7.R1.fq.gz
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN5_S7.R2.fq.gz
  IgG:
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN6_S8.R1.fq.gz
    - /temp_work/ch228298/CUTRUN_0315/Fastq/CUTandRUN6_S8.R2.fq.gz

3DChromatin_ReplicateQC/install_scripts/install_3DChromatin_ReplicateQC.sh --pathtopython /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/envs/3DChromatin_ReplicateQC/bin/python\
 --pathtor /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/envs/3DChromatin_ReplicateQC/bin/R \
 --rlib /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/envs/3DChromatin_ReplicateQC/lib/R/library 

for cond in KD shNT
do
    for treat in plusCL minusCL
    do
        for rep in 1 2 3 4
        do
            source /programs/biogrids.shrc
            name=${cond}_${treat}_rep${rep}
            echo $name
            mkdir -p $name
            trim_galore --paired -o ${name} Fastq/${name}.R1.fq.gz Fastq/${name}.R2.fq.gz

            prefix=${name}/${name}
            fq1=${name}/${name}.R1_val_1.fq.gz
            fq2=${name}/${name}.R2_val_2.fq.gz

            nThreads=8

            ref='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/bwa_GRCm38/GRCm38.primary_assembly.genome.fa'
            genome_file='/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/GRCm38.primary_assembly.genome.fa.Genome'

            bwa mem -t $nThreads -SP5M $ref $fq1 $fq2 | samtools view -Shb - > ${prefix}.bam

            source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate

            mkdir -p ${name}/temp1
            pairtools parse --add-columns mapq --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $nThreads --nproc-out $nThreads --chroms-path $genome_file ${prefix}.bam | pairtools sort --tmpdir=${name}/temp1/ --nproc 8 |pairtools dedup --nproc-in $nThreads --nproc-out $nThreads --mark-dups --output-stats ${prefix}_stats.txt|pairtools split --nproc-in $nThreads --nproc-out $nThreads --output-pairs ${prefix}.pairs
            echo "+++++Finished"
        done
    done
done

for cond in KD shNT
do
    for treat in plusCL minusCL
    do
        for rep in 1 2 
        do
            source /programs/biogrids.shrc
            name=${cond}_${treat}_rep${rep}
            echo $name
            mv Fastq/${name}.R1.fq.gz Fastq/${name}.R1.fq
            gzip Fastq/${name}.R1.fq
            mv Fastq/${name}.R2.fq.gz Fastq/${name}.R2.fq
            gzip Fastq/${name}.R2.fq
            echo "+++++Finished"
        done
    done
done


pairtools merge -o ${name}.pairs --nproc 16 ../KD_minusCL_rep1/KD_minusCL_rep1.pairs ../KD_minusCL_rep2/KD_minusCL_rep2.pairs ../KD_minusCL_rep3/KD_minusCL_rep3.pairs ../KD_minusCL_rep4/KD_minusCL_rep4.pairs


python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/mustache/mustache/diff_mustache.py -f1 microc_KD_plusCL_merge_mapq5.hic -f2 microc_shNT_minusCL_merge_mapq5.hic -r 1000 -st 0.88 -o diff_mustache_plusCL_vs_minusCL/diff_mustache_plusCL_vs_minusCL_1k -pt 0.1 -pt2 0.1 -p 8
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/mustache/mustache/diff_mustache.py -f1 microc_KD_plusCL_merge_mapq5.hic -f2 microc_shNT_minusCL_merge_mapq5.hic -r 5000 -st 0.88 -o diff_mustache_plusCL_vs_minusCL/diff_mustache_plusCL_vs_minusCL_5k -pt 0.1 -pt2 0.1 -p 8
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/mustache/mustache/diff_mustache.py -f1 microc_KD_plusCL_merge_mapq5.hic -f2 microc_shNT_minusCL_merge_mapq5.hic -r 10000 -st 0.88 -o diff_mustache_plusCL_vs_minusCL/diff_mustache_plusCL_vs_minusCL_10k -pt 0.1 -pt2 0.1 -p 8
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/mustache/mustache/diff_mustache.py -f1 microc_KD_plusCL_merge_mapq5.hic -f2 microc_shNT_minusCL_merge_mapq5.hic -r 25000 -st 0.88 -o diff_mustache_plusCL_vs_minusCL/diff_mustache_plusCL_vs_minusCL_25k -pt 0.1 -pt2 0.1 -p 8
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/mustache/mustache/diff_mustache.py -f1 microc_shNT_plusCL_merge_mapq5.hic -f2 microc_shNT_minusCL_merge_mapq5.hic -r 50000 -st 0.88 -o diff_mustache_plusCL_vs_minusCL/diff_mustache_plusCL_vs_minusCL_50k -pt 0.1 -pt2 0.1 -p 8
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/mustache/mustache/diff_mustache.py -f1 microc_shNT_plusCL_merge_mapq5.hic -f2 microc_shNT_minusCL_merge_mapq5.hic -r 100000 -st 0.88 -o diff_mustache_plusCL_vs_minusCL/diff_mustache_plusCL_vs_minusCL_100k -pt 0.1 -pt2 0.1 -p 8

meta$new_celltype = meta$cell_type
meta[meta[,'cell_type'] == 'Pdgfra_APC','new_celltype'] = meta[meta[,'cell_type'] == 'Pdgfra_APC','cell_type_cluster']



source("https://bioconductor.org/biocLite.R")
#biocLite("hicrep", lib=.libPaths()[1]) #old hicrep
biocLite("rhdf5")

install.packages("rmarkdown",lib=.libPaths()[1],repos="http://cran.rstudio.com/")
install.packages("testthat",lib=.libPaths()[1],repos="http://cran.rstudio.com/")
install.packages("reshape2",lib=.libPaths()[1],repos="http://cran.rstudio.com/")
install.packages("pheatmap",lib=.libPaths()[1],repos="http://cran.rstudio.com/")

# hicrep from the hicrep paper
# install.packages("Supplemental_hicrep_1.0.1.tar.gz",dependencies="logical")

#install newest hicrep
biocLite("hicrep", lib=.libPaths()[1])


for i in KD shNT
do 
    for j in minusCL plusCL
    do
        mv ${i}_${j}_ATAC.rep1_peaks.narrowPeak ${i}_${j}_rep3.rep1_peaks.narrowPeak
    done
done


 wiq_result/H2AZ_vs_IgG_pooled_adipo_H2A.Z_ChIP_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig
 wiq_result/H2AZ_CL_vs_IgG_pooled_thermo_adipo_H2A.Z_ChIP_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig


pairtools merge -o ${pair_new} --nproc 16 ${name}_rep1.pairs.mapq5.pairs \
${name}_rep2.pairs.mapq5.pairs \
${name}_rep3.pairs.mapq5.pairs \
${name}_rep4.pairs.mapq5.pairs \
${name}_rep1_pilot.pairs.mapq5.pairs \
${name}_rep2_pilot.pairs.mapq5.pairs \
${name}_rep3_pilot.pairs.mapq5.pairs \
${name}_rep4_pilot.pairs.mapq5.pairs 

source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate
for i in shNT KD
do
    for j in minusCL plusCL
    do
        for n in 1 2 3 4
        do
            name=${i}_${j}_rep${n} #'shNT_plusCL_rep1'
            echo ++$name
            pair_new=${name}/${name}.pairs.mapq5.pairs
            chromsize=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize

            echo "+++convert to cool"
            for r in 2500000 1000000 500000 250000 100000 50000 25000 10000 5000 2000 1000
            do
                echo ++++${r}
                cooler_path=${pair_new}_${r}.cool
                cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 ${chromsize}:${r} ${pair_new} ${cooler_path}
                echo ++++balance
                cooler balance ${cooler_path}
            done
        done
    done
done



import os,sys
import pandas as pd

path=sys.argv[1]

files = [os.path.join(path, x) for x in os.listdir(path) if x.endswith('_50000.txt.gz')]
for f in files:
    d = pd.read_csv(f, header = None, sep = '\t', compression = 'gzip')
    d[0] = d[0].str.rstrip('chr')
    d[2] = d[2].str.rstrip('chr')
    d.to_csv(f+'.new.gz', header = None, index = None, sep = '\t', compression = 'gzip')


for i in `ls *.txt.gz`; do zcat $i |sed 's/chr//g' |gzip > ${i}.new.gz; done
for i in `ls *.txt.gz`; do mv -f ${i}.new.gz ${i}; done



python cellphonedb_run.py CCI_datasets/human_paad/PDAC_A_ST1/sc_meta.tsv CCI_datasets/human_paad/PDAC_A_ST1/sc_count.tsv cellphonedb_output_human_PDAC_A_ST1
python cellphonedb_run.py CCI_datasets/human_paad/PDAC_B_ST3/sc_meta.tsv CCI_datasets/human_paad/PDAC_B_ST3/sc_counts.tsv cellphonedb_output_human_PDAC_B_ST3
python cellphonedb_run.py CCI_datasets/human_scc/P10_rep1_GSM4284325/sc_meta.tsv CCI_datasets/human_scc/P10_rep1_GSM4284325/sc_counts.tsv cellphonedb_output_human_scc_P10_rep1_GSM4284325
python cellphonedb_run.py CCI_datasets/human_scc/P2_rep2_GSM4284317/sc_meta.tsv CCI_datasets/human_scc/P2_rep2_GSM4284317/sc_count.tsv cellphonedb_output_human_scc_P2_rep2_GSM4284317
python cellphonedb_run.py CCI_datasets/human_scc/P5_rep3_GSM4284321/sc_meta.tsv CCI_datasets/human_scc/P5_rep3_GSM4284321/sc_counts.tsv cellphonedb_output_human_scc_P5_rep3_GSM4284321



echo "++++ activating"
source /programs/biogrids.shrc


computeMatrix reference-point -S H2AZ_vs_IgG_pooled_adipo_H2A.Z_ChIP_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig.bw H2AZ_CL_vs_IgG_pooled_thermo_adipo_H2A.Z_ChIP_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig.bw \
                             -R UpTss_shNT_plusCL_vs_minusCL.bed -a 5000 -b 5000 \
                             --outFileName ChIP_2021_UpTss_shNT_plusCL_vs_minusCL_Tss5kb.mat.gz --numberOfProcessors 16

plotHeatmap -m ChIP_2021_UpTss_shNT_plusCL_vs_minusCL_Tss5kb.mat.gz \
      -out ChIP_2021_UpTss_shNT_plusCL_vs_minusCL_Tss5kb.mat.png \
      --colorMap Reds \
      --refPointLabel 'center' --xAxisLabel 'distance (bp)' --yAxisLabel 'peaks' \
      --samplesLabel minusCL plusCL --kmeans 3



computeMatrix reference-point -S H2AZ_vs_IgG_pooled_adipo_H2A.Z_ChIP_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig.bw H2AZ_CL_vs_IgG_pooled_thermo_adipo_H2A.Z_ChIP_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig.bw \
                             -R DownTss_shNT_plusCL_vs_minusCL.bed -a 5000 -b 5000 \
                             --outFileName ChIP_2021_DownTss_shNT_plusCL_vs_minusCL_Tss5kb.mat.gz --numberOfProcessors 16

plotHeatmap -m ChIP_2021_DownTss_shNT_plusCL_vs_minusCL_Tss5kb.mat.gz \
      -out ChIP_2021_DownTss_shNT_plusCL_vs_minusCL_Tss5kb.mat.png \
      --colorMap Reds \
      --refPointLabel 'center' --xAxisLabel 'distance (bp)' --yAxisLabel 'peaks' \
      --samplesLabel minusCL plusCL --kmeans 3
echo "++++ Finished"


import os,sys
import _bw
import numpy as np
import pandas as pd

d=[1000, 10000, 30000, 100000]
path='/temp_work/ch228298/H2AZ_Danpos/qnor_bsub_bw'
bwFiles = os.listdir(path)

tss = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.protein_coding.tss.csv'
rp_res = {}
for f in bwFiles:
    label=f.split('.sorted.bam')[0]
    f = os.path.join(path, f)
    for dd in d:
        if str(dd) not in rp_res:
            rp_res[str(dd)] = {}
        _bw.getrp(f, tss, label+'_lisaRP_%sDecay.txt'%dd, dd, 0, 0)
        rp_res[label]=pd.read_csv(label+'_lisaRP_%sDecay.txt'%dd, sep = '\t', header = None)

for dd in d:
    rp_res_mat = pd.DataFrame(map(lambda x: rp_res[dd][x][4], rp_res[dd]))
    rp_res_mat = rp_res_mat.T
    rp_res_mat.index = rp_res[label][3].tolist()
    rp_res_mat.columns = rp_res.keys()
    rp_res_mat.to_csv('combined_lisaRP_orignal_%sDecay.csv'%dd)


d1 = read.csv('secretions.tsv', sep = '\t', row.names = 1)
d2 = read.csv('../../mBAT_cold2/secretions.tsv', sep = '\t', row.names = 1)


cor.test(d1[grepl('e]$', rownames(d1)),1], d2[grepl('e]$', rownames(d2)),1])



awk '{
    if ( $0 ~ /^#chromsize/ ) {
        if ( $0 ~ /chr[0-9]|chrX|chrY|chrM/ ) {
            print
        }
    }else if ( $0 ~ /^#samheader/ && $0 ~ /SQ/) {
        if ( $0 ~ /chr[0-9]|chrX|chrY|chrM/ ) {
            print
        }
    }else if ( $0 ~ /^#/ ){
        print
    }else{
        if ( $2 ~ /chr[0-9]|chrX|chrY|chrM/ && $4 ~ /chr[0-9]|chrX|chrY|chrM/ && $9 > 5 && $10 > 5 && (($3 - $5) > 100 || ($3 - $5) < -100)){
            print
        }
    }
}' ${pair} > test.pair

awk 'if ( $2 ~ /chr[0-9]|chrX|chrY|chrM/ && $4 ~ /chr[0-9]|chrX|chrY|chrM/ && $9 > 5 && $10 > 5 && $(($3 >= $5 ? $3 - $5 : $5 - $3)) < 100){
            print
        }' test.pair > test2.pair

awk '{print $(($3 >= $5 ? $3 - $5 : $5 - $3))}' test2.pair 
awk '{print $($3 - $5) }' test2.pair 

  KD_plusCL_H2AZ_abcam:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H2AZ-abcam_S12_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H2AZ-abcam_S12_R2.fq.gz
  KD_plusCL_H2AZ_AM_rep1:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H2AZ-Active-Motif-rep1_S2_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H2AZ-Active-Motif-rep1_S2_R2.fq.gz
  KD_plusCL_H2AZ_AM_rep2:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H2AZ-Active-Motif-rep2_S10_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H2AZ-Active-Motif-rep2_S10_R2.fq.gz
  KD_plusCL_H3K27ac_rep1:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K27ac-rep1_S4_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K27ac-rep1_S4_R2.fq.gz
  KD_plusCL_H3K27ac_rep2:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K27ac-rep2_S14_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K27ac-rep2_S14_R2.fq.gz
  KD_plusCL_H3K4me3_rep1:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K4me3-rep1_S6_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K4me3-rep1_S6_R2.fq.gz
  KD_plusCL_H3K4me3_rep2:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K4me3-rep2_S16_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-H3K4me3-rep2_S16_R2.fq.gz
  KD_plusCL_IgG:
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-IgG_S8_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/KD-CL-IgG_S8_R2.fq.gz
  shNT_plusCL_H2AZ_abcam:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H2AZ-abcam_S11_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H2AZ-abcam_S11_R2.fq.gz
  shNT_plusCL_H2AZ_AM_rep1:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H2AZ-Active-Motif-rep1_S1_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H2AZ-Active-Motif-rep1_S1_R2.fq.gz
  shNT_plusCL_H2AZ_AM_rep2:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H2AZ-Active-Motif-rep2_S9_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H2AZ-Active-Motif-rep2_S9_R2.fq.gz
  shNT_plusCL_H3K27ac_rep1:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K27ac-rep1_S3_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K27ac-rep1_S3_R2.fq.gz
  shNT_plusCL_H3K27ac_rep2:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K27ac-rep2_S13_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K27ac-rep2_S13_R2.fq.gz
  shNT_plusCL_H3K4me3_rep1:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K4me3-rep1_S5_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K4me3-rep1_S5_R2.fq.gz
  shNT_plusCL_H3K4me3_rep2:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K4me3-rep2_S15_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-H3K4me3-rep2_S15_R2.fq.gz
  shNT_plusCL_IgG:
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-IgG_S7_R1.fq.gz
    - /temp_work/ch228298/ChIP_202306/ready/NT-CL-IgG_S7_R2.fq.gz



sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/Danpos_bash.sh ../ChIP_202306/analysis/align/KD_plusCL_H2AZ_abcam/KD_plusCL_H2AZ_abcam_unique.sorted.bam ../ChIP_202306/analysis/align/KD_plusCL_IgG/KD_plusCL_IgG_unique.sorted.bam ChIP_202306_KD_plusCL_H2AZabcam_vs_IgG
sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/Danpos_bash.sh ../ChIP_202306/analysis/align/KD_plusCL_H2AZ_AM_rep1/KD_plusCL_H2AZ_AM_rep1_unique.sorted.bam ../ChIP_202306/analysis/align/KD_plusCL_IgG/KD_plusCL_IgG_unique.sorted.bam ChIP_202306_KD_plusCL_H2AZ_AM_vs_IgG_rep1
sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/Danpos_bash.sh ../ChIP_202306/analysis/align/KD_plusCL_H2AZ_AM_rep2/KD_plusCL_H2AZ_AM_rep2_unique.sorted.bam ../ChIP_202306/analysis/align/KD_plusCL_IgG/KD_plusCL_IgG_unique.sorted.bam ChIP_202306_KD_plusCL_H2AZ_AM_vs_IgG_rep2
sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/Danpos_bash.sh ../ChIP_202306/analysis/align/shNT_plusCL_H2AZ_abcam/shNT_plusCL_H2AZ_abcam_unique.sorted.bam ../ChIP_202306/analysis/align/shNT_plusCL_IgG/shNT_plusCL_IgG_unique.sorted.bam ChIP_202306_shNT_plusCL_H2AZabcam_vs_IgG
sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/Danpos_bash.sh ../ChIP_202306/analysis/align/shNT_plusCL_H2AZ_AM_rep1/shNT_plusCL_H2AZ_AM_rep1_unique.sorted.bam ../ChIP_202306/analysis/align/shNT_plusCL_IgG/shNT_plusCL_IgG_unique.sorted.bam ChIP_202306_shNT_plusCL_H2AZ_AM_vs_IgG_rep1
sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/CommonCode/Danpos_bash.sh ../ChIP_202306/analysis/align/shNT_plusCL_H2AZ_AM_rep2/shNT_plusCL_H2AZ_AM_rep2_unique.sorted.bam ../ChIP_202306/analysis/align/shNT_plusCL_IgG/shNT_plusCL_IgG_unique.sorted.bam ChIP_202306_shNT_plusCL_H2AZ_AM_vs_IgG_rep2

### plot heatmap

source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate
chromsize=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize

for pair in `ls ../../../*_rep*/*.pairs.mapq5.pairs`
do
    echo "+++filter mapq"
    pair_new=`basename ${pair}`
    awk '{
        if ( $0 ~ /^#chromsize/ ) {
            if ( $0 ~ /chr[0-9]|chrX|chrY|chrM/ ) {
                print
            }
        }else if ( $0 ~ /^#samheader/ && $0 ~ /SQ/) {
            if ( $0 ~ /chr[0-9]|chrX|chrY|chrM/ ) {
                print
            }
        }else if ( $0 ~ /^#/ ){
            print
        }else{
            if ( $2 ~ /chr[0-9]|chrX|chrY|chrM/ && $4 ~ /chr[0-9]|chrX|chrY|chrM/ && $9 > 5 && $10 > 5 && (($3 - $5) > 100 || ($3 - $5) < -100)){
                print
            }
        }
    }' ${pair} > ${pair_new}

    echo "+++convert to cool"
    for r in 2500000 1000000 500000 250000 100000 50000 25000 10000 5000 2000 1000
    do
        echo ++++${r}
        cooler_path=cool/${pair_new}_${r}.cool
        cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 ${chromsize}:${r} ${pair_new} ${cooler_path}
        echo ++++balance
        cooler balance ${cooler_path}
    done

d1 = pd.read_csv('shNT_plusCL_mapq5_merge.pairs_5000.bedpe', sep = '\t', header = None)
d2 = pd.read_csv('shNT_minusCL_mapq5_merge.pairs_5000.bedpe', sep = '\t', header = None)
d1_loop = d1[d1[6] > 0.98]
d2_loop = d2[d2[6] > 0.98]

d1_index = d1_loop[0]+':'+d1_loop[1].astype('str')+':'+d1_loop[2].astype('str')+':'+d1_loop[3]+':'+d1_loop[4].astype('str')+':'+d1_loop[5].astype('str')
d2_index = d2_loop[0]+':'+d2_loop[1].astype('str')+':'+d2_loop[2].astype('str')+':'+d2_loop[3]+':'+d2_loop[4].astype('str')+':'+d2_loop[5].astype('str')

d = pd.merge(d1[[6,7]], d2[[6,7]], left_index = True, right_index = True)


for i in adipose bone_marrow breast_milk liver pancreas adrenal_gland brain kidney lung skin
do
echo '#!/bin/bash' > ${i}_comp.sbatch
echo "#SBATCH --partition=cbp-compute # queue to be used" >> ${i}_comp.sbatch
echo "#SBATCH --time=100:00:00 # Running time (in hours-minutes-seconds)" >> ${i}_comp.sbatch
echo "#SBATCH --job-name=${i} # Job name" >> ${i}_comp.sbatch
echo "#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails" >> ${i}_comp.sbatch
echo "#SBATCH --mail-user=your_email_address # Email address to send the job status" >> ${i}_comp.sbatch
echo "#SBATCH --output=${i}.txt # Name of the output file" >> ${i}_comp.sbatch
echo "#SBATCH --nodes=1 # Number of compute nodes" >> ${i}_comp.sbatch
echo "#SBATCH --ntasks=12 # Number of cpu cores on one node" >> ${i}_comp.sbatch
echo "#SBATCH --mem=50G" >> ${i}_comp.sbatch
echo "source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate" >> ${i}_comp.sbatch
echo "module load singularity" >> ${i}_comp.sbatch
echo "echo 'runing compass'" >> ${i}_comp.sbatch
t=${i}
exp_tsv=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/DISCO/${t}_avg_exp.tsv
out_dir=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/DISCO/compass_res/${t}_res
temp_dir=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/DISCO/compass_res/${t}_temp
echo "singularity run --bind /lab-share/Cardio-Chen-e2/Public/rongbinzheng:/lab-share/Cardio-Chen-e2/Public/rongbinzheng /lab-share/Cardio-Chen-e2/Public/rongbinzheng/tmp/compass sh /lab-share/Cardio-Chen-e2/Public/rongbinzheng/metabolism/test_compass/compass_run.sh ${exp_tsv} homo_sapiens ${out_dir} ${temp_dir} 12" >> ${i}_comp.sbatch
echo "echo 'Finished'" >> ${i}_comp.sbatch
done


https://zenodo.org/record/7381371/files/disco_adipose_v01.h5ad?download=1
https://zenodo.org/record/7381371/files/disco_bone_marrow_v01.h5ad?download=1


compass_met_ann = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/Compass/compass/Resources/Recon2_export/met_md.csv')

compass_rxn_ann = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/Compass/compass/Resources/Recon2_export/rxn_md.csv')

met_ann = pd.read_csv('/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/MEBOCOST/data/mebocost_db/common/metabolite_annotation_HMDB_summary.tsv',
                     sep = '\t')
alias = {str(i).upper():[str(x).upper() for x in str(j).split('; ')] for i,j in met_ann[['metabolite', 'synonyms_name']].values.tolist()}


python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/peakachu/diffPeakachu/pair-probs.py shNT_plusCL_mapq5_merge.pairs_5000_loops.0.95.bedpe shNT_minusCL_mapq5_merge.pairs_5000_loops.0.95.bedpe shNT_plusCL-minusCL_union.loops.0.95.bedpe
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/peakachu/diffPeakachu/diffPeakachu.py shNT_plusCL_mapq5_merge.pairs_5000_loops.0.95.bedpe shNT_minusCL_mapq5_merge.pairs_5000_loops.0.95.bedpe shNT_plusCL-minusCL_union.loops.0.95.bedpe

python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/peakachu/diffPeakachu/pair-probs.py shNT_plusCL_mapq5_merge.pairs_5000_loops.0.98.bedpe shNT_minusCL_mapq5_merge.pairs_5000_loops.0.98.bedpe shNT_plusCL-minusCL_union.loops.0.98.bedpe
python /lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/peakachu/diffPeakachu/diffPeakachu.py shNT_plusCL_mapq5_merge.pairs_5000_loops.0.98.bedpe shNT_minusCL_mapq5_merge.pairs_5000_loops.0.98.bedpe shNT_plusCL-minusCL_union.loops.0.98.bedpe

wig1=ChIP_04082022_H2AZ_plusCL_unique_rep1.sorted.bam.dedup.bgsub.Fnor.wig 
wig2=../wiq_result/ChIP_04082022_H2AZ_minusCL_unique_rep1.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
label=ChIP_04082022_plusCL_vs_minusCL_H2AZ_rep1
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dregion ${wig1}:${wig2} -o qnorm/dregion/${label}
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${wig1}:${wig2} -o qnorm/dpeak/${label}

wig1=../wiq_result/ChIP_04082022_H2AZ_plusCL_unique_rep2.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
wig2=../wiq_result/ChIP_04082022_H2AZ_minusCL_unique_rep2.sorted.bam.dedup.bgsub.Fnor.qnor.wig
label=ChIP_04082022_plusCL_vs_minusCL_H2AZ_rep2
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dregion ${wig1}:${wig2} -o qnorm/dregion/${label}
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${wig1}:${wig2} -o qnorm/dpeak/${label}

wig1=../wiq_result/ChIP_04082022_H2AZ_plusCL_unique_rep3.sorted.bam.dedup.bgsub.Fnor.qnor.wig  
wig2=../wiq_result/ChIP_04082022_H2AZ_minusCL_unique_rep3.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
label=ChIP_04082022_plusCL_vs_minusCL_H2AZ_rep3
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dregion ${wig1}:${wig2} -o qnorm/dregion/${label}
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${wig1}:${wig2} -o qnorm/dpeak/${label}

wig1=../wiq_result/ChIP_04282023_H2AZ_ab4174_plusCL_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
wig2=../wiq_result/ChIP_04282023_H2AZ_ab4174_minusCL_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
label=ChIP_04282023_H2AZ_ab4174_plusCL_vs_minusCL_H2AZ_rep1
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dregion ${wig1}:${wig2} -o qnorm/dregion/${label}
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${wig1}:${wig2} -o qnorm/dpeak/${label}

wig1=../wiq_result/ChIP_04282023_H2AZ_AM_plusCL_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
wig2=../wiq_result/ChIP_04282023_H2AZ_AM_minusCL_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
label=ChIP_04282023_H2AZ_AM_plusCL_vs_minusCL_H2AZ_rep1
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dregion ${wig1}:${wig2} -o qnorm/dregion/${label}
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${wig1}:${wig2} -o qnorm/dpeak/${label}

wig1=CUTRUN_01162023_H2AZ_plusCL_rep1_unique.sorted.bam.dedup.bgsub.Fnor.wig 
wig2=../wiq_result/CUTRUN_01162023_H2AZ_minusCL_rep1_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
label=CUTRUN_01162023_H2AZ_ab4174_plusCL_vs_minusCL_H2AZ_rep1
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dregion ${wig1}:${wig2} -o qnorm/dregion/${label}
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${wig1}:${wig2} -o qnorm/dpeak/${label}

wig1=../wiq_result/CUTRUN_01162023_H2AZ_plusCL_rep2_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
wig2=../wiq_result/CUTRUN_01162023_H2AZ_minusCL_rep2_unique.sorted.bam.dedup.bgsub.Fnor.qnor.wig 
label=CUTRUN_01162023_H2AZ_ab4174_plusCL_vs_minusCL_H2AZ_rep2
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dregion ${wig1}:${wig2} -o qnorm/dregion/${label}
python.danpos2 /programs/x86_64-linux/danpos2/2.2.2/danpos.py dpeak ${wig1}:${wig2} -o qnorm/dpeak/${label}


Rscript neurochat_run.R CCI_datasets/human_heart/sc_normed.tsv CCI_datasets/human_heart/heart_sc_meta.tsv human_heart
Rscript neurochat_run.R CCI_datasets/mouse_cortex/sc_normed.tsv CCI_datasets/mouse_cortex/sc_meta.tsv mouse_cortex
Rscript neurochat_run.R CCI_datasets/human_intestinal/ST_A3_GSM4797918/sc_normed.tsv CCI_datasets/human_intestinal/ST_A3_GSM4797918/sc_meta.tsv human_intestinal_A3
Rscript neurochat_run.R CCI_datasets/human_intestinal/ST_A4_GSM4797919/sc_normed.tsv CCI_datasets/human_intestinal/ST_A4_GSM4797919/sc_meta.tsv human_intestinal_A4
Rscript neurochat_run.R CCI_datasets/human_paad/PDAC_A_ST1/sc_normed.tsv CCI_datasets/human_paad/PDAC_A_ST1/sc_meta.tsv human_PDAC_A
Rscript neurochat_run.R CCI_datasets/human_paad/PDAC_B_ST3/sc_normed.tsv CCI_datasets/human_paad/PDAC_B_ST3/sc_meta.tsv human_PDAC_B
Rscript neurochat_run.R CCI_datasets/human_scc/P10_rep1_GSM4284325/sc_normed.tsv CCI_datasets/human_scc/P10_rep1_GSM4284325/sc_meta.tsv human_SCC_P10
Rscript neurochat_run.R CCI_datasets/human_scc/P2_rep2_GSM4284317/sc_normed.tsv CCI_datasets/human_scc/P2_rep2_GSM4284317/sc_meta_TN_Epithelial.tsv human_SCC_P2
Rscript neurochat_run.R CCI_datasets/human_scc/P5_rep3_GSM4284321/sc_normed.tsv CCI_datasets/human_scc/P5_rep3_GSM4284321/sc_meta.tsv human_SCC_P5



for i in `ls DU145_AAVS1_dmso*`; do sbatch -A cbp $i; done
for i in `ls DU145_AAVS1_VE822*`; do sbatch -A cbp $i; done
for i in `ls DU145_KPNA3_dmso*`; do sbatch -A cbp $i; done
for i in `ls DU145_KPNA3_VE822*`; do sbatch -A cbp $i; done
for i in `ls 22RV1_AAVS1_dmso*`; do sbatch -A cbp $i; done
for i in `ls 22RV1_AAVS1_VE822*`; do sbatch -A cbp $i; done
for i in `ls 22RV1_KPNA3_dmso*`; do sbatch -A cbp $i; done
for i in `ls 22RV1_KPNA3_VE822*`; do sbatch -A cbp $i; done

#!/bin/bash
#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=5:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=run2 # Job name
#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address # Email address to send the job status
#SBATCH --output=output_%A_%a.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=16 # Number of cpu cores on one node
#SBATCH --mem=20G

echo "++++ activating"
source /programs/biogrids.shrc

bed=
computeMatrix reference-point -S /temp_work/ch228298/ChIP_202306/analysis/peaks/shNT_plusCL_H3K27ac.rep1/shNT_plusCL_H3K27ac.rep1_treat_pileup.bw \
                            /temp_work/ch228298/ChIP_202306/analysis/peaks/shNT_plusCL_H3K27ac.rep2/shNT_plusCL_H3K27ac.rep2_treat_pileup.bw \
                            /temp_work/ch228298/ChIP_202306/analysis/peaks/shNT_plusCL_H3K4me3.rep1/shNT_plusCL_H3K4me3.rep1_treat_pileup.bw \
                            /temp_work/ch228298/ChIP_202306/analysis/peaks/shNT_plusCL_H3K4me3.rep1/shNT_plusCL_H3K4me3.rep1_treat_pileup.bw \
                             -R $bed -a 2000 -b 2000 \
                             --outFileName ${bed}_K27_K4.mat.gz --numberOfProcessors 16


samples:
  minusCL_FSK_input:
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/1_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/1_R2_001.fastq.gz
  minusCL_FSK_H2AZ:  
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/2_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/2_R2_001.fastq.gz
  minusCL_FSK_H3K27ac:  
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/3_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/3_R2_001.fastq.gz
  minus_FSK_H3K4me3:  
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/4_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/4_R2_001.fastq.gz
  plusCL_FSK_input:  
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/5_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/5_R2_001.fastq.gz
  plusCL_FSK_H2AZ:  
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/6_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/6_R2_001.fastq.gz
  plusCL_FSK_H3K27ac:  
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/7_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/7_R2_001.fastq.gz
  plusCL_FSK_H3K4me3:  
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/8_R1_001.fastq.gz
    - /lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/H2AZ_human/30-962665206/00_fastq/8_R2_001.fastq.gz
    
minusCL_FSK_H2A.Z,minusCL_FSK_H2A.Z,minusCL_FSK_input,,
minusCL_FSK_H3K27ac,minusCL_FSK_H3K27ac,minusCL_FSK_input,,
minusCL_FSK_H3K4me3,minusCL_FSK_H3K4me3,minusCL_FSK_input,,
plusCL_FSK_H2A.Z,plusCL_FSK_H2A.Z,plusCL_FSK_input,,
plusCL_FSK_H3K27ac,plusCL_FSK_H3K27ac,plusCL_FSK_input,,
plusCL_FSK_H3K4me3,plusCL_FSK_H3K4me3,plusCL_FSK_input,,


r=5000
for i in 1 2 3 4
do
    rclone copy -P google-drive:"TsengLab/Shared/H2AZ_BrownAdipo/MicroC_Analysis_2023/shNT_plusCL_rep"${i}"/shNT_plusCL_rep"${i}".pairs.mapq5.pairs_"${r}".cool" .
    rclone copy -P google-drive:"TsengLab/Shared/H2AZ_BrownAdipo/MicroC_Analysis_2023/shNT_minusCL_rep"${i}"/shNT_minusCL_rep"${i}".pairs.mapq5.pairs_"${r}".cool" .
    rclone copy -P google-drive:"TsengLab/Shared/H2AZ_BrownAdipo/MicroC_Analysis_2023/KD_plusCL_rep"${i}"/KD_plusCL_rep"${i}".pairs.mapq5.pairs_"${r}".cool" .
done



bamCoverage --outFileName plusCL_H2AZ_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20220408_H2AZ_ChIP_Yang/analysis/H2AZ_2022_1st/analysis/align/H2AZ_CL/H2AZ_CL_unique.sorted.bam
bamCoverage --outFileName plusCL_IgG_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20220408_H2AZ_ChIP_Yang/analysis/H2AZ_2022_1st/analysis/align/IgG_CL/IgG_CL_unique.sorted.bam
bamCoverage --outFileName plusCL_H2AZ_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20220408_H2AZ_ChIP_Yang/analysis/H2AZ_2022_2st/analysis/align/H2AZ_CL/H2AZ_CL_unique.sorted.bam
bamCoverage --outFileName plusCL_IgG_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20220408_H2AZ_ChIP_Yang/analysis/H2AZ_2022_2st/analysis/align/IgG_CL/IgG_CL_unique.sorted.bam
bamCoverage --outFileName plusCL_H2AZ_rep3.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20220408_H2AZ_ChIP_Yang/analysis/H2AZ_2022_3st/analysis/align/H2AZ_CL/H2AZ_CL_unique.sorted.bam
bamCoverage --outFileName plusCL_IgG_rep3.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20220408_H2AZ_ChIP_Yang/analysis/H2AZ_2022_3st/analysis/align/IgG_CL/IgG_CL_unique.sorted.bam


bamCoverage --outFileName minusCL_H3K27ac.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/H3K27ac_minusCL/H3K27ac_minusCL_unique.sorted.bam
bamCoverage --outFileName plusCL_H3K27ac.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/H3K27ac_plusCL/H3K27ac_plusCL_unique.sorted.bam
bamCoverage --outFileName minusCL_H3K4me3.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/H3K4me3_minusCL/H3K4me3_minusCL_unique.sorted.bam
bamCoverage --outFileName plusCL_H3K4me3.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/H3K4me3_plusCL/H3K4me3_plusCL_unique.sorted.bam

bamCoverage --outFileName minusCL_IgG.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/IgG_minusCL/IgG_minusCL_unique.sorted.bam
bamCoverage --outFileName plusCL_IgG.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/IgG_plusCL/IgG_plusCL_unique.sorted.bam
bamCoverage --outFileName minusCL_Input.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/Input_minusCL/Input_minusCL_unique.sorted.bam
bamCoverage --outFileName plusCL_Input.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230428_ChIPseq_Day11/processed/analysis/align/Input_plusCL/Input_plusCL_unique.sorted.bam

bamCoverage --outFileName shNT_plusCL_H3K27ac_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/shNT_plusCL_H3K27ac_rep1/shNT_plusCL_H3K27ac_rep1_unique.sorted.bam
bamCoverage --outFileName shNT_plusCL_H3K27ac_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/shNT_plusCL_H3K27ac_rep2/shNT_plusCL_H3K27ac_rep2_unique.sorted.bam
bamCoverage --outFileName shNT_plusCL_H3K4me3_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/shNT_plusCL_H3K4me3_rep1/shNT_plusCL_H3K4me3_rep1_unique.sorted.bam
bamCoverage --outFileName shNT_plusCL_H3K4me3_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/shNT_plusCL_H3K4me3_rep2/shNT_plusCL_H3K4me3_rep2_unique.sorted.bam
bamCoverage --outFileName KD_plusCL_H3K27ac_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/KD_plusCL_H3K27ac_rep1/KD_plusCL_H3K27ac_rep1_unique.sorted.bam
bamCoverage --outFileName KD_plusCL_H3K27ac_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/KD_plusCL_H3K27ac_rep2/KD_plusCL_H3K27ac_rep2_unique.sorted.bam
bamCoverage --outFileName KD_plusCL_H3K4me3_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/KD_plusCL_H3K4me3_rep1/KD_plusCL_H3K4me3_rep1_unique.sorted.bam
bamCoverage --outFileName KD_plusCL_H3K4me3_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/KD_plusCL_H3K4me3_rep2/KD_plusCL_H3K4me3_rep2_unique.sorted.bam
bamCoverage --outFileName KD_plusCL_IgG.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/KD_plusCL_IgG/KD_plusCL_IgG_unique.sorted.bam
bamCoverage --outFileName shNT_plusCL_IgG.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230616_ChIP_H2AZ_Histone/processed/analysis/align/shNT_plusCL_IgG/shNT_plusCL_IgG_unique.sorted.bam


bamCoverage --outFileName KD_plusCL_ATAC_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20221209_MicroC_ATAC_Yang/20221209_ATAC1rep/analysis/align/KD_plusCL_ATAC/KD_plusCL_ATAC_unique.sorted.bam
bamCoverage --outFileName shNT_plusCL_ATAC_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20221209_MicroC_ATAC_Yang/20221209_ATAC1rep/analysis/align/shNT_plusCL_ATAC/shNT_plusCL_ATAC_unique.sorted.bam
bamCoverage --outFileName shNT_minusCL_ATAC_rep2.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20221209_MicroC_ATAC_Yang/20221209_ATAC1rep/analysis/align/shNT_minusCL_ATAC/shNT_minusCL_ATAC_unique.sorted.bam

bamCoverage --outFileName KD_plusCL_ATAC_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230106_0107_ATACseq_Fastq_regenerated/processed/analysis/align/KD_plusCL_rep2/KD_plusCL_rep2_unique.sorted.bam
bamCoverage --outFileName shNT_minusCL_ATAC_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230106_0107_ATACseq_Fastq_regenerated/processed/analysis/align/shNT_minusCL_rep2/shNT_minusCL_rep2_unique.sorted.bam
bamCoverage --outFileName shNT_plusCL_ATAC_rep1.bw --outFileFormat bigwig --normalizeUsing CPM --numberOfProcessors 12 --bam 20230106_0107_ATACseq_Fastq_regenerated/processed/analysis/align/shNT_plusCL_rep2/shNT_plusCL_rep2_unique.sorted.bam

mv KD_plusCL_rep1.R1.fq.gz KD_plusCL_MicroC_rep1.R1.fq.gz
mv KD_plusCL_rep1.R2.fq.gz KD_plusCL_MicroC_rep1.R2.fq.gz
mv KD_plusCL_rep2.R1.fq.gz KD_plusCL_MicroC_rep2.R1.fq.gz
mv KD_plusCL_rep2.R2.fq.gz KD_plusCL_MicroC_rep2.R2.fq.gz
mv KD_plusCL_rep3.R1.fq.gz KD_plusCL_MicroC_rep3.R1.fq.gz
mv KD_plusCL_rep3.R2.fq.gz KD_plusCL_MicroC_rep3.R2.fq.gz
mv KD_plusCL_rep4.R1.fq.gz KD_plusCL_MicroC_rep4.R1.fq.gz
mv KD_plusCL_rep4.R2.fq.gz KD_plusCL_MicroC_rep4.R2.fq.gz

mv MicroC-10_S3.R1.fq.gz shNT_plusCL_MicroC_rep3.R1.fq.gz
mv MicroC-10_S3.R2.fq.gz shNT_plusCL_MicroC_rep3.R2.fq.gz
mv MicroC-14_S2.R1.fq.gz shNT_plusCL_MicroC_rep4.R1.fq.gz
mv MicroC-14_S2.R2.fq.gz shNT_plusCL_MicroC_rep4.R2.fq.gz
mv MicroC-2_S1.R1.fq.gz shNT_plusCL_MicroC_rep1.R1.fq.gz
mv MicroC-2_S1.R2.fq.gz shNT_plusCL_MicroC_rep1.R2.fq.gz
mv MicroC-6_S2.R1.fq.gz shNT_plusCL_MicroC_rep2.R1.fq.gz
mv MicroC-6_S2.R2.fq.gz shNT_plusCL_MicroC_rep2.R2.fq.gz
mv MicroC-13_S1.R1.fq.gz shNT_minusCL_MicroC_rep4.R1.fq.gz
mv MicroC-13_S1.R2.fq.gz shNT_minusCL_MicroC_rep4.R2.fq.gz
mv MicroC-9_S1.R1.fq.gz shNT_minusCL_MicroC_rep3.R1.fq.gz
mv MicroC-9_S1.R2.fq.gz shNT_minusCL_MicroC_rep3.R2.fq.gz
mv shNT-minusCL-rep1_S1.R1.fq.gz shNT_minusCL_MicroC_rep1.R1.fq.gz
mv shNT-minusCL-rep1_S1.R2.fq.gz shNT_minusCL_MicroC_rep1.R2.fq.gz
mv shNT-minusCL-rep2_S2.R1.fq.gz shNT_minusCL_MicroC_rep2.R1.fq.gz
mv shNT-minusCL-rep2_S2.R2.fq.gz shNT_minusCL_MicroC_rep2.R2.fq.gz


mv KD_plusCL_rep3.R1.fq.gz KD_plusCL_MicroC_pilot_rep3.R1.fq.gz
mv KD_plusCL_rep3.R2.fq.gz KD_plusCL_MicroC_pilot_rep3.R2.fq.gz
mv KD_plusCL_rep4.R1.fq.gz KD_plusCL_MicroC_pilot_rep4.R1.fq.gz
mv KD_plusCL_rep4.R2.fq.gz KD_plusCL_MicroC_pilot_rep4.R2.fq.gz
mv shNT_minusCL_rep3.R1.fq.gz shNT_minusCL_MicroC_pilot_rep3.R1.fq.gz
mv shNT_minusCL_rep3.R2.fq.gz shNT_minusCL_MicroC_pilot_rep3.R2.fq.gz
mv shNT_minusCL_rep4.R1.fq.gz shNT_minusCL_MicroC_pilot_rep4.R1.fq.gz
mv shNT_minusCL_rep4.R2.fq.gz shNT_minusCL_MicroC_pilot_rep4.R2.fq.gz
mv shNT_plusCL_rep3.R1.fq.gz shNT_plusCL_MicroC_pilot_rep3.R1.fq.gz
mv shNT_plusCL_rep3.R2.fq.gz shNT_plusCL_MicroC_pilot_rep3.R2.fq.gz
mv shNT_plusCL_rep4.R1.fq.gz shNT_plusCL_MicroC_pilot_rep4.R1.fq.gz
mv shNT_plusCL_rep4.R2.fq.gz shNT_plusCL_MicroC_pilot_rep4.R2.fq.gz


mv KD-CL-H3K27ac-rep1_S4_R1.fq.gz KD_plusCL_H3K27ac_rep1.R1.fq.gz 
mv KD-CL-H3K27ac-rep1_S4_R2.fq.gz KD_plusCL_H3K27ac_rep1.R2.fq.gz 
mv KD-CL-H3K27ac-rep2_S14_R1.fq.gz KD_plusCL_H3K27ac_rep2.R1.fq.gz 
mv KD-CL-H3K27ac-rep2_S14_R2.fq.gz KD_plusCL_H3K27ac_rep2.R2.fq.gz 
mv KD-CL-H3K4me3-rep1_S6_R1.fq.gz KD_plusCL_H3K4me3_rep1.R1.fq.gz 
mv KD-CL-H3K4me3-rep1_S6_R2.fq.gz KD_plusCL_H3K4me3_rep1.R2.fq.gz 
mv KD-CL-H3K4me3-rep2_S16_R1.fq.gz KD_plusCL_H3K4me3_rep2.R1.fq.gz 
mv KD-CL-H3K4me3-rep2_S16_R2.fq.gz KD_plusCL_H3K4me3_rep2.R2.fq.gz 
mv KD-CL-IgG_S8_R1.fq.gz KD_plusCL_IgG.R1.fq.gz 
mv KD-CL-IgG_S8_R2.fq.gz KD_plusCL_IgG.R2.fq.gz 
mv NT-CL-H3K27ac-rep1_S3_R1.fq.gz shNT_plusCL_H3K27ac_rep1.R1.fq.gz 
mv NT-CL-H3K27ac-rep1_S3_R2.fq.gz shNT_plusCL_H3K27ac_rep1.R2.fq.gz 
mv NT-CL-H3K27ac-rep2_S13_R1.fq.gz shNT_plusCL_H3K27ac_rep2.R1.fq.gz 
mv NT-CL-H3K27ac-rep2_S13_R2.fq.gz shNT_plusCL_H3K27ac_rep2.R2.fq.gz 
mv NT-CL-H3K4me3-rep1_S5_R1.fq.gz shNT_plusCL_H3K4me3_rep1.R1.fq.gz 
mv NT-CL-H3K4me3-rep1_S5_R2.fq.gz shNT_plusCL_H3K4me3_rep1.R2.fq.gz 
mv NT-CL-H3K4me3-rep2_S15_R1.fq.gz shNT_plusCL_H3K4me3_rep2.R1.fq.gz 
mv NT-CL-H3K4me3-rep2_S15_R2.fq.gz shNT_plusCL_H3K4me3_rep2.R2.fq.gz 
mv NT-CL-IgG_S7_R1.fq.gz shNT_plusCL_IgG.R1.fq.gz 
mv NT-CL-IgG_S7_R2.fq.gz shNT_plusCL_IgG.R2.fq.gz 


mv ChIP_04082022_H2AZ_minusCL_unique_rep1.sorted.bam.dedup.bgsub.Fnor.peaks.xls minusCL_H2AZ_rep1_peaks.xls
mv ChIP_04082022_H2AZ_minusCL_unique_rep2.sorted.bam.dedup.bgsub.Fnor.peaks.xls minusCL_H2AZ_rep2_peaks.xls
mv ChIP_04082022_H2AZ_minusCL_unique_rep3.sorted.bam.dedup.bgsub.Fnor.peaks.xls minusCL_H2AZ_rep3_peaks.xls
mv ChIP_04082022_H2AZ_plusCL_unique_rep1.sorted.bam.dedup.bgsub.Fnor.peaks.xls plusCL_H2AZ_rep1_peaks.xls
mv ChIP_04082022_H2AZ_plusCL_unique_rep2.sorted.bam.dedup.bgsub.Fnor.peaks.xls plusCL_H2AZ_rep2_peaks.xls
mv ChIP_04082022_H2AZ_plusCL_unique_rep3.sorted.bam.dedup.bgsub.Fnor.peaks.xls plusCL_H2AZ_rep3_peaks.xls


for i in 1 2 3 5 6
do
rclone copyto -P google-drive:Collaboration/Burns/RNAseq/KO7_rep${i}/analysis/STAR_outputAligned.sortedByCoord.out.bam.bai KO7_rep${i}.bam.bai
done

for i in 2 5 6
do
rclone copyto -P google-drive:Collaboration/Burns/RNAseq/WT_rep${i}/analysis/STAR_outputAligned.sortedByCoord.out.bam.bai WT_rep${i}.bam.bai
done


for i in 1 2 3 4 6
do
rclone copyto -P google-drive:Collaboration/Burns/RNAseq/KO21_rep${i}/analysis/STAR_outputAligned.sortedByCoord.out.bam KO21_rep${i}.bam
done


rmatsshi=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/rmats2sashimiplot/src/rmats2sashimiplot/rmats2sashimiplot.py

python ${rmatsshi} --b1 WT.txt --b2 KO7.txt --event-type MXE -e ./plot_events/RBFOX2_KO7_vs_WT_MXE.txt --l1 WT --l2 KO7 --exon_s 1 --intron_s 5 -o sashimiplot_MXE_group --group-info grouping.gf
python ${rmatsshi} --b1 WT.txt --b2 KO7.txt --event-type MXE -e ./plot_events/RBFOX2_KO7_vs_WT_MXE.txt --l1 WT --l2 KO7 --exon_s 1 --intron_s 5 -o sashimiplot_MXE --fig-height 12

python ${rmatsshi} --b1 WT.txt --b2 KO7.txt --event-type SE -e ./plot_events/RBFOX2_KO7_vs_WT_SE.txt --l1 WT --l2 KO7 --exon_s 1 --intron_s 5 -o sashimiplot_SE_group --group-info grouping.gf
python ${rmatsshi} --b1 WT.txt --b2 KO7.txt --event-type SE -e ./plot_events/RBFOX2_KO7_vs_WT_SE.txt --l1 WT --l2 KO7 --exon_s 1 --intron_s 5 -o sashimiplot_SE --fig-height 12

fq1=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/00_fastq/1_R1_001.fastq.gz
fq2=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/00_fastq/1_R2_001.fastq.gz
ad=ACGGACTT
n=NT_minusCL_bisox
out=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/trimed_fq
trim_galore --paired --dont_gzip -o $out --basename $n --fastqc_args '-d ${QC}' -j 8 --adapter $ad $fq1 $fq2

cutadapt --version
4.8

n=NT_minusCL_bisox
fq1=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/trimed_fq/${n}_val_1.fq
fq2=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/trimed_fq/${n}_val_2.fq
fa=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/ref_files/mm10/bwa_indices/mm10.fa
m10size=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize
DISPbind align -i $fa -n $n -a $fq1 -b $fq2 -o processed_out -p 8 -g $m10size

nm=NT_plusCL
bam1=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_DMSO.sorted.bam
bam2=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_bisox.sorted.bam
macs2 callpeak --nomodel -B --SPMR -f BAMPE -g mm -t $bam2 -c $bam1 --outdir callpeak/ -n ${nm}_bisox_peaks --broad --keep-dup 1

nm=NT_minusCL
bam1=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_DMSO.sorted.bam
bam2=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_bisox.sorted.bam
macs2 callpeak --nomodel -B --SPMR -f BAMPE -g mm -t $bam2 -c $bam1 --outdir callpeak/ -n ${nm}_bisox_peaks --broad --keep-dup 1

nm=KD_minusCL
bam1=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_DMSO.sorted.bam
bam2=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_bisox.sorted.bam
macs2 callpeak --nomodel -B --SPMR -f BAMPE -g mm -t $bam2 -c $bam1 --outdir callpeak/ -n ${nm}_bisox_peaks --broad --keep-dup 1

nm=KD_plusCL
bam1=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_DMSO.sorted.bam
bam2=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/processed_out/${nm}_bisox.sorted.bam
macs2 callpeak --nomodel -B --SPMR -f BAMPE -g mm -t $bam2 -c $bam1 --outdir callpeak/ -n ${nm}_bisox_peaks --broad --keep-dup 1


for n in NT_minusCL_DMSO NT_plusCL_bisox NT_plusCL_DMSO KD_minusCL_bisox KD_minusCL_DMSO KD_plusCL_bisox KD_plusCL_DMSO
do
    echo '#!/bin/bash' > ${n}_map.sbatch
    echo '#SBATCH --partition=cbp-compute # queue to be used' >> ${n}_map.sbatch
    echo '#SBATCH --time=50:00:00 # Running time (in hours-minutes-seconds)' >> ${n}_map.sbatch
    echo '#SBATCH --job-name=disp # Job name' >> ${n}_map.sbatch
    echo '#SBATCH --mail-type=BEGIN,END,FAIL # send and email when the job begins, ends or fails' >> ${n}_map.sbatch
    echo '#SBATCH --mail-user=your_email_address # Email address to send the job status' >> ${n}_map.sbatch
    echo '#SBATCH --output=output_%A_%a.txt # Name of the output file' >> ${n}_map.sbatch
    echo '#SBATCH --nodes=1 # Number of compute nodes' >> ${n}_map.sbatch
    echo '#SBATCH --ntasks=8 # Number of cpu cores on one node' >> ${n}_map.sbatch
    echo '#SBATCH --mem=30G' >> ${n}_map.sbatch

    echo "echo '++++ activating'" >> ${n}_map.sbatch
    echo "source /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/activate" >> ${n}_map.sbatch
    echo "conda activate chips" >> ${n}_map.sbatch
    echo "fq1=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/trimed_fq/${n}_val_1.fq" >> ${n}_map.sbatch
    echo "fq2=/temp_work/ch228298/Nature_revise_Yang/DisP_seq_062024/trimed_fq/${n}_val_2.fq" >> ${n}_map.sbatch
    echo "fa=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/software/ref_files/mm10/bwa_indices/mm10.fa" >> ${n}_map.sbatch
    echo "m10size=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize" >> ${n}_map.sbatch
    cmd="DISPbind align -i $fa -n $n -a $fq1 -b $fq2 -o processed_out -p 8 -g $m10size"
    echo $cmd >> ${n}_map.sbatch
    echo 'echo "Finished"' >> ${n}_map.sbatch
done


6 genes have insufficient target regions and may have compromised detection efficiency: 
TNP1
TUBA3C
AP002518.2
CPT1B
CHKB
PAOX
1 gene is untargetable: TPSB2


1. setup.py vs. requirements.txt
2. https://stackoverflow.com/questions/75682385/runtimeerror-cuda-error-no-kernel-image-is-available-for-execution-on-the-devi

#slice_ids = slice_ids[::-1]

args = easydict.EasyDict({})
args.epochs = 40
args.lr = 0.001
args.k = 20
args.alpha = 0.1 # weight of transcriptional loss
args.diff_omics = False # whether to use different omics data
args.mode = 'None' # Choose the mode among 'align', 'stitch' and None
args.dimension = 2  # choose the dimension of coordinates (2 or 3)
args.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

src = sc.read_h5ad('src.h5ad')
tgt = sc.read_h5ad('tgt.h5ad')
sc.pp.subsample(src, 0.1)
sc.pp.subsample(tgt, 0.1)

### align adjacent slices by SANTO
aligned_src_cor, trans_dict = santo(src,tgt, args)
cor = np.array(src.obsm['spatial'])
cor_tgt = np.array(tgt.obsm['spatial'])
cor_new = np.dot(cor, trans_dict['fine_R_ab'].T) + trans_dict['fine_T_ab']
pcc,ci = evaluation(cor_new,
                     cor_tgt,
                     src.X,
                     tgt.X,
                     np.array(src.obs.cell_label_updated),
                     np.array(tgt.obs.cell_label_updated))
print(f'accuracy: {pcc}, ci: {ci}')


res = NULL
for (r in c('250000', '100000', '50000', '10000')){
    tmp = list()
    for (i in paste0('chr', c(as.character(1:19)))){
        d = readRDS(file = paste0('hicrep_reslu_', r, '_', i, '.rds'))
        tmp[[i]] = do.call(c, lapply(d, function(x) {x[['scc']]}))
    }
    tmp = as.data.frame(do.call(rbind, tmp))
    tmp$chrom = rownames(tmp)
    tmp$resolu = rep(r, nrow(tmp))
    res = rbind(res, tmp)
}
write.table(res, file = 'replicates_hicrep_corr_mat.txt', sep = '\t', quote = F)

chromsize=/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/mm10.chromSize
for wig in `ls *CTCF_vs_input/pooled/*bgsub.Fnor.wig`
do
    /lab-share/Cardio-Chen-e2/Public/rongbinzheng/anaconda3/bin/wigToBigWig -clip ${wig} $chromsize ${wig}.bw
done

