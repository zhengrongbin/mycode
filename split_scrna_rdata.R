library(Seurat)

d = readRDS('seurat_collapsed_clusters.Rdata')

conds = c('RT', 'TN', 'cold2', 'cold7')

print('==== rds')
for (cond in conds) {
	print(cond)
	tmp = subset(seurat_bulk, `sample` == cond)
	saveRDS(tmp, file = paste0(cond, '_seurat.rds'))
}

print('==== exp')
for (cond in conds) {
	print(cond)
	tmp = subset(seurat_bulk, `sample` == cond)
	exp_mat = as.matrix(tmp@assays$RNA@data)
	write.csv(exp_mat, file = paste0(cond, '_mat.csv'))
}
