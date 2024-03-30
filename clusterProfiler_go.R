args = commandArgs(T)
library(clusterProfiler)
library(ggplot2)
#library(org.Hs.eg.db)
library(org.Mm.eg.db)

input=args[1] # the input file, which contains 1 column gene symbol
label=args[2] # the name label for this running
background = args[3] # either all, or a file with 1 column gene symbol

geneset = as.vector(read.csv(input, header = F, sep = ',')[,1])
geneset = geneset[!grepl('^MT-', geneset)]

GO = list()
for (k in c('BP', 'MF', 'CC')){
  print(k)
  if (background == 'all'){
    GO[[k]] = enrichGO(gene = geneset, OrgDb = org.Mm.eg.db, pvalueCutoff = 1, qvalueCutoff = 1, ont = k, keyType = 'SYMBOL')
  }else{
    background_genes = as.vector(read.csv(background, header = F, sep = '\t')[,1])
    GO[[k]] = enrichGO(gene = geneset, OrgDb = org.Mm.eg.db, pvalueCutoff = 1, qvalueCutoff = 1, ont = k, universe = background_genes, keyType = 'SYMBOL')
  }
}

### KEGG
#geneid = as.vector(bitr(geneset, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Mm.eg.db, drop = TRUE)[,2])
#k='KEGG'
#if (background == 'all'){
#    GO[[k]] = enrichKEGG(gene = geneid, organism='mmu', keyType='kegg', pvalueCutoff = 1, qvalueCutoff = 1)
#  }else{
#    background_genes = as.vector(read.csv(background, header = F, sep = '\t')[,1])
#    GO[[k]] = enrichKEGG(gene = geneid, organism='mmu', keyType='kegg', pvalueCutoff = 1, qvalueCutoff = 1, universe = background_genes)
# }
#
pdf(paste0(label, '.pdf'), width = 10, height = 5)
clusterProfiler::dotplot(GO[['MF']], showCategory = 20, title = paste0(label, '_MF'))
clusterProfiler::dotplot(GO[['BP']], showCategory = 20, title = paste0(label, '_BP'))
clusterProfiler::dotplot(GO[['CC']], showCategory = 20, title = paste0(label, '_CC'))
#clusterProfiler::dotplot(GO[['KEGG']], showCategory = 20, title = paste0(label, '_KEGG'))
dev.off()
saveRDS(GO, paste0(label, '_go.rds'))

