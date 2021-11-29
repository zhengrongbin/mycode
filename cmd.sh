docker pull ubuntu:latest   

## load whole document
docker run -it --name routine -p 8888:8888 -v /Users/rongbinzheng/Documents:/Users/rongbinzheng/Documents routine:latest

## for MECOM
docker run -it --name MECOM -p 8880:8880 -v /Users/rongbinzheng/Documents/BCH/ChenLab:/home/ChenLab -v /Users/rongbinzheng/Documents/CommonData:/home/CommonData seurat:v2.0.1


jupyter notebook --allow-root --port=8888 --no-browser --ip=0.0.0.0

remotes::install_version(package = 'Seurat', version = package_version('2.0.1'))


mkdir /lab-share/Cardio-Chen-e2/Public/rongbinzheng/google_drive
rclone mount google-drive:~/ /project/RC_Cardio-Chen-e2/rongbinzheng/google_drive --allow-other --allow-non-empty --vfs-cache-mode writes


srun -A bch -p bch-interactive --pty /bin/bash


### ========== MGHPCC ========== 
## interactive
srun -A bch-mghpcc --pty /bin/bash


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


















