docker pull ubuntu:latest   

curl -H "Accept: application/vnd.github+json" \
 -H "Authorization: zhengrongbin <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>" \
 https://api.github.com/repos/zhengrongbin/MEBOCOST/traffic/clones


## load whole document
docker run -it --name routine -p 8080:8080 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents routine:latest

docker run -it --name seurat -p 8080:8080 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents seurat:v3.1.5

## for MECOM
docker run -it --name MECOM -p 8880:8880 -v /Users/rongbinzheng/Documents/BCH/ChenLab:/home/ChenLab -v /Users/rongbinzheng/Documents/CommonData:/home/CommonData seurat:v2.0.1

## for compass
docker run -it --name compass -p 8880:8880 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents  compass:latest

## for vue
docker run -it --name vue -p 8880:8880 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents  vue:latest

jupyter notebook --allow-root --port=8080 --no-browser --ip=0.0.0.0 

remotes::install_version(package = 'Seurat', version = package_version('2.0.1'))

## E2
jupyter notebook --no-browser --port 8900 --NotebookApp.iopub_data_rate_limit=10000000000 --ip=0.0.0.0 --NotebookApp.allow_origin=* --allow-root

ssh -N -L 8900:compute-1-6.tch.harvard.edu:8900 -o ServerAliveInterval=30 ch228298@e2.tch.harvard.edu

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

#hicrep from the hicrep paper
#install.packages("Supplemental_hicrep_1.0.1.tar.gz",dependencies="logical")

#install newest hicrep
biocLite("hicrep", lib=.libPaths()[1])


for i in KD shNT
do 
    for j in minusCL plusCL
    do
        mv ${i}_${j}_ATAC_unique.sorted.bam ${i}_${j}_rep3_unique.sorted.bam
    done
done




