library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

# test
# gene = 'KMT2D'
# cancer = 'STAD'
# filepath = paste0('/Users/zhengrongbin/Documents/Workplace/Projects/TAD/diffexp/top_mutated/', 'no_', cancer, '_mut_', gene, '_diffCond1VsCond2.csv')
# diff_exp = read.csv(filepath, row.names = 1)
# diff_exp$tad_label = protein_annotation[as.vector(rownames(diff_exp)),'tad']
# diff_exp_good = diff_exp[!is.na(diff_exp$tad_label),]
# diff_exp_good$tad_cluster = TAD_cluster[as.vector(diff_exp_good$tad_label),'tad_cluster']
# diff_exp_good$geneType = protein_annotation[as.vector(rownames(diff_exp_good)),'gene_label']
# diff_exp_good_coding = subset(diff_exp_good, geneType == 'protein_coding')
# up = subset(diff_exp_good_coding, logFC >= log2(1.5) & adj.P.Val <= 0.05)
# up = up[order(-up$logFC),]
# up_c5_genes = subset(up, tad_cluster == 'cluster_5')
# down = subset(diff_exp_good_coding, logFC <= (-log2(1.5)) & adj.P.Val <= 0.05)
# down = down[order(down$logFC),]
# down_c5_genes = subset(down, tad_cluster == 'cluster_5')

## do GO
GO = function(geneList, geneType = 'SYMBOL', background = FALSE){
  res = list()
  for (t in c('MF', 'BP', 'CC')){
    ego <- enrichGO(gene          = geneList,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = geneType,
                    ont           = t,
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 1)
    res[t] = ego
  }
}

geneList = up_c5_genes$logFC; names(geneList) = rownames(up_c5_genes)
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              keyType = 'SYMBOL',
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              verbose      = FALSE)




