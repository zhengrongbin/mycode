import scanpy as sc
import commot as ct
import pandas as pd

adata = sc.read_h5ad('/temp_work/ch228298/merfish/BAT_Combined_3cold_2TN.h5ad')
adata.obs['cellid'] = [x.split('-')[0] for x in adata.obs.index.tolist()]
adata.obs['celllabel'] = adata.obs['cellid'].astype('str')+'~'+adata.obs['Datasets'].astype('str').str.replace('_region_0', '')
adata.obs['Datasets'] = adata.obs['Datasets'].str.replace('_region_0', '')

## add spatial cordinates
spatial_cord1 = pd.read_csv('/temp_work/ch228298/merfish/output/202307211201_ColdBAT71_VMSC02901/cell_metadata.csv', index_col = 0)
spatial_cord1['celllabel'] = spatial_cord1.index.astype('str')+'~'+'202307211201_ColdBAT71_VMSC02901'

spatial_cord2 = pd.read_csv('/temp_work/ch228298/merfish/output/202307241349_BATTN1_VMSC02901/cell_metadata.csv', index_col = 0)
spatial_cord2['celllabel'] = spatial_cord2.index.astype('str')+'~'+'202307241349_BATTN1_VMSC02901'

spatial_cord3 = pd.read_csv('/temp_work/ch228298/merfish/output/202311011242_coldbatsample2_VMSC02901/region_0//cell_metadata.csv', index_col = 0)
spatial_cord3['celllabel'] = spatial_cord3.index.astype('str')+'~'+'202311011242_coldbatsample2_VMSC02901'

spatial_cord4 = pd.read_csv('/temp_work/ch228298/merfish/output/202401081212_BATcold3_VMSC02901/region_0//cell_metadata.csv', index_col = 0)
spatial_cord4['celllabel'] = spatial_cord4.index.astype('str')+'~'+'202401081212_BATcold3_VMSC02901'

spatial_cord5 = pd.read_csv('/temp_work/ch228298/merfish/output/202401091250_TN2BAT_VMSC02901/region_0//cell_metadata.csv', index_col = 0)
spatial_cord5['celllabel'] = spatial_cord5.index.astype('str')+'~'+'202401091250_TN2BAT_VMSC02901'

spatial_cord = pd.concat([spatial_cord1, spatial_cord2, spatial_cord3, spatial_cord4, spatial_cord5])
spatial_cord.index = spatial_cord['celllabel'].tolist()

## give spatial
adata.obsm['spatial'] = spatial_cord.loc[adata.obs['celllabel'].tolist()][['center_x','center_y','min_x', 'min_y','max_x','max_y']].to_numpy()

df_cellchat = ct.pp.ligand_receptor_database(species='mouse', database='CellChat', signaling_type=None)

### sample by sample
for x in adata.obs['Datasets'].unique().tolist():
	adata1 = adata[adata.obs['Datasets'] == x]
	df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata1, min_cell_pct=0.05)
	ct.tl.spatial_communication(adata1,
	    database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
	adata1.write('/temp_work/ch228298/merfish/commt_res/'+x+'.commt.h5ad')


