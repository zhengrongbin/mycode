args=commandArgs(T)
library(data.table)

exp_file = args[1]
meta_file = args[2]
label = args[3]

exp_mat = fread(exp_file, sep = '\t')
genes = as.vector(exp_mat$V1)
exp_mat = as.data.frame(exp_mat[,2:ncol(exp_mat)])
rownames(exp_mat) = genes
meta_mat = as.data.frame(read.table(meta_file, sep = '\t', row.names = 1))

inter_cell = intersect(rownames(meta_mat), colnames(exp_mat))
exp_mat = exp_mat[,inter_cell]
meta_mat = meta_mat[inter_cell,,drop=F]

library(NeuronChat)

x <- createNeuronChat(as.matrix(exp_mat), DB='human', group.by = meta_mat[,1]);
x <- run_NeuronChat(x,M=1000)
net_aggregated_x <- net_aggregation(x@net,method = 'weight')

comm_res = NULL
for (n in names(x@net0)){
    tmp = as.data.frame(reshape2::melt(x@net0[[n]]))
    colnames(tmp) = c('sender', 'receiver', 'strength')
    tmp$lr = n
    comm_res = rbind(comm_res, tmp)
}
comm_res = as.data.frame(comm_res)

pvalue = NULL
for (n in 1:length(x@pvalue)){
    tmp = as.data.frame(reshape2::melt(x@pvalue[[n]]))
    pvalue = c(pvalue, tmp$value)
}
comm_res$pval = pvalue

db = do.call(rbind, x@DB)
comm_res_new = merge(comm_res, db, by.x = 'lr', by.y = 0)

write.table(comm_res_new, file = paste0('NeuroChat_res/', label, '_neurochat.tsv'), 
    sep = '\t', quote = F, row.names = F)
