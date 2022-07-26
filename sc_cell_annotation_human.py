import pandas as pd
import numpy as np
import scanpy as sc

## ================
## scbat_obj_paga is the scanpy object
## this example I am doing human single cell RNA-seq data cell annotation
## ==============


## ==== prepare maker genes from PanglaoDB
pangodb = pd.read_csv('PanglaoDB_markers_27_Mar_2020.tsv', sep = '\t')
# homology = pd.read_csv('human_mouse_homology_gene_pair.csv')
# ## to mouse
# marker_gene_mat = pd.merge(homology, marker_gene_mat, left_on = 'human_gene', right_index = True)
# marker_gene_mat.index = marker_gene_mat['mouse_gene'].tolist()
# marker_gene_mat = marker_gene_mat.drop(['mouse_gene', 'human_gene'], axis = 1)
### select species, here for mouse
marker_gene_mat = pangodb[(pangodb['species'].isin(['Hs', 'Mm Hs']))]
# ## only need some tissue and cell types
# tissue_need = ['Connective tissue', 'Vasculature', 'Immune system', 'Smooth muscle']
# celltype_need = ['Schwann cells', 'Neurons']
# marker_gene_mat = marker_gene_mat.loc[(marker_gene_mat['organ'].isin(tissue_need)) |
#                                       (marker_gene_mat['cell type'].isin(celltype_need)) ,:]
# ## remove very rare cell types that may cause misunderstanding in annotation
del_cells = ['Chondrocytes', 'Gamma delta T cells', 'Plasmacytoid dendritic cells', 'Red pulp macrophages',
            'Endothelial cells (aorta)', 'Myofibroblasts', 'Myoepithelial cells']
marker_gene_mat = marker_gene_mat[~marker_gene_mat['cell type'].isin(del_cells)]

## === select genes with high sensitivity (frequency of gene expressed in the cell type) and low specificity (frequency of gene expressed in not of the cells) score
scbat_obj_paga_genes = scbat_obj_paga_update.raw.to_adata().var.index.tolist()

marker_genes_dict = {}
for celltype in marker_gene_mat['cell type'].unique().tolist():
    tmp = marker_gene_mat[marker_gene_mat['cell type'] == celltype]
    tmp.index = tmp['official gene symbol'].tolist()
#     tmp = tmp['sensitivity_mouse']
#     tmp = tmp[['sensitivity_mouse', 'specificity_mouse']].T.mean().dropna().sort_values()
    ttmp = tmp[(tmp['sensitivity_human'] > 0.1) & (tmp['specificity_human'] < 0.4)].index.tolist()
#     ttmp = homology[homology['human_gene'].isin(tmp.index.tolist())]['mouse_gene'].tolist()
    # ttmp = homology[homology['human_gene'].isin(ttmp)]['mouse_gene'].tolist() ## since the gene name in Pangolao DB always upper characters, and gene can be homology
    ttmp = list(set(ttmp) & set(scbat_obj_paga_genes))
    if ttmp:
        marker_genes_dict[celltype] = ttmp
    
  
## update marker genes manually
marker_genes_dict['Vascular smooth muscle cells'].extend(['Acta2', 'Pdgfrb', 'Cspg4'])
marker_genes_dict['ASCs'] = ['Cd34', 'Pdgfra', 'Ly6a', 'Itgb1']
marker_genes_dict['Preadipocytes'] = ['Cd34', 'Pdgfra', 'Ly6a', 'Pparg', 'Cebpa', 'Dcn']
marker_genes_dict['Adipocytes'].extend(['Retn', 'Cidec', 'Ucp1'])

## update Cd4 and Cd8 T cell
marker_genes_dict['CD8 T cells'] = marker_genes_dict['T cells'] + ['Cd8a', 'Cd8b1']
marker_genes_dict['CD4 T cells'] = marker_genes_dict['T cells'] + ['Cd4']
marker_genes_dict.pop('T cells')
### T reg
marker_genes_dict['T regulatory cells'].extend(['Foxp3', 'Il2ra'])

## correct: Nuocyte has been called as ILC2s since 2013, https://en.wikipedia.org/wiki/ILC2#cite_note-pmid24876829-9
## and add two additional markers
marker_genes_dict['ILC2s'] = marker_genes_dict['Nuocytes'] + ['Il5', 'Il13']
marker_genes_dict.pop('Nuocytes')

## remove cell type with only one gene as marker
marker_count = {x:len(marker_genes_dict[x]) for x in marker_genes_dict}
for cell in [x for x in marker_count if marker_count[x] < 5]:
    marker_genes_dict.pop(cell)
# remove T memory, since get same marker from T cells after checking, keep only T cell
marker_genes_dict.pop('T memory cells')
    
### visulize by bubal plot
sc.pl.dotplot(scbat_obj_paga, marker_genes_dict, 'cluster.new', dendrogram=True,
              save = 'cell_marker_dot.pdf')

## extract diff score for all genes and make matrix
cluster_diff_score_mat = pd.DataFrame()
for c in scbat_obj_paga.obs['cluster.new'].unique().tolist():
    tmp = sc.get.rank_genes_groups_df(scbat_obj_paga, c)
    tmp.columns = c + '_' + tmp.columns
    tmp.index = tmp[c+'_names'].tolist()
    cluster_diff_score_mat = pd.concat([cluster_diff_score_mat, tmp[[c+'_scores']]], axis = 1)
cluster_diff_score_mat.columns = cluster_diff_score_mat.columns.str.replace('_scores', '')

## marker gene score 
cell_cluster_marker_score = {}
for cell_type in marker_genes_dict:
#     print(cell_type)
    genes = marker_genes_dict[cell_type]
    cell_cluster_marker_score[cell_type] = cluster_diff_score_mat.loc[genes].mean()
cell_cluster_marker_score = pd.DataFrame(cell_cluster_marker_score).T

df1 = cell_cluster_marker_score.copy()
df1.index = [x+' (n=%s)'%marker_count[x] for x in df1.index.tolist()]
df1.columns.name = 'Cluster'
## this is adipose tissue data, so remove Hepatic stellate cells (n=13) and Pancreatic stellate cells (n=15)
## this usually do by manual
df1 = df1.drop(['Hepatic stellate cells (n=13)', 'Pancreatic stellate cells (n=15)'])

### select best annotation by the score and plot heatmap
best_selected_celltype = df1.apply(lambda x: x.index[x==x.max()]).T[0].unique().tolist()
df1_new = df1.reindex(best_selected_celltype)

sp = sns.clustermap(df1_new,
            cmap = 'bwr', figsize = (16, 8), linewidth = .5, center = 0,
              vmax = 50, vmin = -50, annot = df1_new, fmt = '.2f',
              annot_kws={"size": 4})
plt.close()

col_order = [x.get_text() for x in sp.ax_heatmap.get_xticklabels()]
row_order = [x.get_text() for x in sp.ax_heatmap.get_yticklabels()]

df1_new = df1_new.reindex(index = row_order, columns = col_order)
fig, ax = plt.subplots(figsize = (18, 8))
sns.heatmap(data = df1_new,
           center = 0,vmax = 50, vmin = -50, annot = df1_new, fmt = '.2f',
              annot_kws={"size": 6}, cmap = 'bwr', cbar_kws = {'shrink':.5},
           linewidth = .1)
ax.set_title('Marker Score', fontsize = 14, pad = 10)
plt.tight_layout()
plt.savefig('cell_marker_score_heatmap_stat.pdf')
plt.show()
plt.close()







