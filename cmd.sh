docker pull ubuntu:latest   

curl -H "Accept: application/vnd.github+json" \
 -H "Authorization: zhengrongbin <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>" \
 https://api.github.com/repos/zhengrongbin/MEBOCOST/traffic/clones


## load whole document
docker run -it --name routine -p 8080:8080 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents routine:latest

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

ssh -N -L 8900:compute-1-0.tch.harvard.edu:8900 -o ServerAliveInterval=30 ch228298@e2.tch.harvard.edu

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


