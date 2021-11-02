#"""GO analysis"""
args = commandArgs(T)

f = args[1] # output of DAVID table
save_f = args[2] # hhe file name you want to save

#f = paste0('./tmp/blood/', 'intersection.txt')
go_table = read.csv(f, header = T, sep = '\t')

go_entry = go_table[which(sapply(as.vector(go_table[,1]), function(x){return(length(grep('GOTERM', x)) > 0)})),]
go_entry = cbind(go_entry, 'adj.p.value'=p.adjust(as.vector(go_entry[,'PValue']), method = 'BH'))
go_entry = go_entry[order(go_entry[,'adj.p.value']),]
kegg_entry = go_table[which(sapply(as.vector(go_table[,1]), function(x){return(length(grep('KEGG', x)) > 0)})),]
kegg_entry = cbind(kegg_entry, 'adj.p.value'=p.adjust(as.vector(kegg_entry[,'PValue']), method = 'BH'))
kegg_entry = kegg_entry[order(kegg_entry[,'adj.p.value']),]
go_entry_top = -log10(go_entry[1:15,'adj.p.value'])
names(go_entry_top) = sapply(as.vector(go_entry[1:15,"Term"]), function(x){return(strsplit(x, '\\~')[[1]][2])})
kegg_entry_top = -log10(kegg_entry[1:15,'adj.p.value'])
names(kegg_entry_top) = sapply(as.vector(kegg_entry[1:15,"Term"]), function(x){return(strsplit(x, '\\:')[[1]][2])})
#pdf('./tmp/blood/intersection.pdf', width=10, height=6)
pdf(save_f, width=10, height=8)
attach(mtcars)
opar= par(no.readonly=TRUE)
par(mfrow=c(2,1), mar = c(6,40,4,4), cex = 0.7)
barplot(go_entry_top, las = 2, horiz = T, col = 'orange', xlab = "-log10(adj.p.value)", main = 'GO', font = 2)
barplot(kegg_entry_top, las = 2, horiz = T, col = 'orange', xlab = "-log10(adj.p.value)", main = 'KEGG', font = 2)
par(opar)
detach(mtcars)
dev.off()

# plot(1:5, 1:5, axes = F, type = 'n', xlab = '', ylab = '')
# ven = venn.diagram(promoter.variance.gene.top[8:10],filename =  NULL, col = c('red', 'blue', 'orange'))
# grid.draw(ven)


# blood = promoter.variance.gene.top[c('B Lymphocyte', 'Monocyte')]
# plot(1:5, 1:5, axes = F, type = 'n', xlab = '', ylab = '')
# ven = venn.diagram(blood, filename =  NULL)
# grid.draw(ven)

# intersect_gene = intersect(blood[[1]], blood[[2]])
# B_cell = setdiff(blood[[1]], blood[[2]])
# Monocyte = setdiff(blood[[2]], blood[[1]])

# for (n in 1:10){
# 	system(paste0('Rscript ../../../../code/GO_code.R 200gene_sampleCluster_',n,'_GO.txt 200gene_sampleCluster_',n,'_GO.pdf'))
# }






