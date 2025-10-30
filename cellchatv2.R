args=base::commandArgs(TRUE)
library(CellChat)
library(data.table)

exp_file = args[1]
meta_file = args[2]
label = args[3]

exp_mat = fread(exp_file, sep = '\t')
genes = as.vector(exp_mat$V1)
exp_mat = as.data.frame(exp_mat[,2:ncol(exp_mat)])
rownames(exp_mat) = genes

meta = read.csv(meta_file, sep = '\t', row.names = 1, header = F)
colnames(meta) <- c('ct')

## normalize count data
# data.input.norm = normalizeData(data.input)
## create cellchat object
cellchat <- createCellChat(object = as.matrix(exp_mat), meta = as.matrix(meta), group.by = "ct")


## set database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
## extract CCC as a data frame
df.net <- subsetCommunication(cellchat)

saveRDS(cellchat, paste0('cellchatv2_res/', label, '.rds'))
write.table(df.net, paste0('cellchatv2_res/', label, '.tsv', sep = '\t'))
