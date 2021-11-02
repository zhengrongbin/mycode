
import os,sys
import _bw
import numpy as np
import pandas as pd

bwFiles = ['../adipo_H2AZ_H2AZac/adipo_H2AZ_H2AZac/peaks/adipo_H2AZac.rep1/adipo_H2AZac.rep1_treat_pileup.bw', 
'../adipo_H2AZ_H2AZac/adipo_H2AZ_H2AZac/peaks/adipo_H2AZ.rep1/adipo_H2AZ.rep1_treat_pileup.bw',
'../thermo_adipo_H2AZ_H2AZac/thermo_adipo_H2AZ_H2AZac/peaks/thermo_adip_H2AZac.rep1/thermo_adip_H2AZac.rep1_treat_pileup.bw',
'../thermo_adipo_H2AZ_H2AZac/thermo_adipo_H2AZ_H2AZac/peaks/thermo_adip_H2AZ.rep1/thermo_adip_H2AZ.rep1_treat_pileup.bw']

tss = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.protein_coding.tss.csv'
rp_res = {}
for f in bwFiles:
	_bw.getrp(f, tss, os.path.basename(f).replace('.rep1_treat_pileup.bw', '_lisaRP_1kDecay.txt'), 1000.0, 0, 0)
	rp_res[os.path.basename(f).replace('.rep1_treat_pileup.bw', '')]=pd.read_csv(os.path.basename(f).replace('.rep1_treat_pileup.bw', '_lisaRP_1kDecay.txt'), sep = '\t', header = None)

rp_res_mat = pd.DataFrame(map(lambda x: rp_res[x][4], rp_res))
rp_res_mat = rp_res_mat.T
rp_res_mat.index = rp_res['adipo_H2AZac'][3].tolist()
rp_res_mat.columns = rp_res.keys()
rp_res_mat.to_csv('combined_lisaRP_orignal_1kDecay.csv')
# ## normalize by sqrt
# rp_res_mat_norm = rp_res_mat.apply(lambda col: np.sqrt(col))

bwFiles = ['/lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/bigwig/adipo_H2AZac.bw',
'/lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/bigwig/adipo_H2AZ.bw',
'/lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/bigwig/thermo_adipo_H2AZac.bw',
'/lab-share/Cardio-Chen-e2/Public/rongbinzheng/DataProcess/bigwig/thermo_adipo_H2AZ.bw']

tss = '/lab-share/Cardio-Chen-e2/Public/rongbinzheng/Genome/mm10/gencode.vM23.annotation.protein_coding.tss.csv'
rp_res = {}
for f in bwFiles:
	_bw.getrp(f, tss, os.path.basename(f).replace('.bw', '_lisaRP_1kDecay.txt'), 1000.0, 0, 0)
	rp_res[os.path.basename(f).replace('.bw', '')]=pd.read_csv(os.path.basename(f).replace('.bw', '_lisaRP_1kDecay.txt'), sep = '\t', header = None)

rp_res_mat = pd.DataFrame(map(lambda x: rp_res[x][4], rp_res))
rp_res_mat = rp_res_mat.T
rp_res_mat.index = rp_res['adipo_H2AZac'][3].tolist()
rp_res_mat.columns = rp_res.keys()
rp_res_mat.to_csv('combined_lisaRP_orignal_1kDecay.csv')
