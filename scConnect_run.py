import os,sys
import pandas as pd
import scanpy as sc
import scConnect as cn
import matplotlib
import matplotlib.pyplot as plt

h5ad_file = sys.argv[1]
label = sys.argv[2]

h5ad_file = 'CCI_datasets/human_heart/scanpy.h5ad'
adata = sc.read_h5ad(h5ad_file)

adata_tissue = cn.genecall.meanExpression(adata, groupby="celltype", normalization=False, use_raw=False, transformation="log1p")
adata_tissue = cn.connect.ligands(adata_tissue)
adata_tissue = cn.connect.receptors(adata_tissue)
adata_tissue = cn.connect.specificity(adata_tissue, n=10, groupby="celltype")
cn.connect.save_specificity(adata_tissue, "%s_celltype_specificity.xlsx"%label)
edges = cn.connect.interactions(emitter=adata_tissue, target=adata_tissue, self_reference=True)
nodes = cn.connect.nodes(adata_tissue)
G_tissue = cn.graph.build_graph(edges, nodes)
Gs_tissue = cn.graph.split_graph(G_tissue)


