git clone https://RongbinZheng@bitbucket.org/plumbers/cidc_chips.git

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

cd cidc_chips  

conda install mamba -n base -c conda-forge

mamba env create -f environment.yml -n chips

conda activate chips

## Homer
perl ~/miniconda3/envs/chips/share/homer/.//configureHomer.pl -install
perl ~/miniconda3/envs/chips/share/homer/.//configureHomer.pl -install hg38
perl ~/miniconda3/envs/chips/share/homer/.//configureHomer.pl -install mm10
