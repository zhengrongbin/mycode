d = read.table('./data/gencode.v19.long_noncoding_RNAs.gtf', sep = '\t')
dgene = subset(d, V3 == 'gene')
dgene_ann = sapply(as.vector(dgene$V9), function(x){
  tmp = do.call(rbind, strsplit(strsplit(as.character(x), '\\; ')[[1]], '\\ '))
  ttmp =as.vector(tmp[,2]); names(ttmp) = as.vector(tmp[,1])
  return(ttmp[c('gene_id', 'transcript_id', 'gene_type', 'gene_name', 'transcript_name')])
})
dgene_ann = t(dgene_ann)
rownames(dgene_ann) = seq(1, nrow(dgene_ann))
dgene_ann = cbind(dgene[,c(1,4,5,7)], dgene_ann)
dgene_ann$length = abs(dgene_ann$V5-dgene_ann$V4)
rownames(dgene_ann) = as.vector(dgene_ann$gene_id)
colnames(dgene_ann) = c('chr', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'gene_type', 'gene_name', 'transcript_name', 'length')
write.table(dgene_ann, './data/gencode19release_lncRNAgene_feature.txt', quote = F, row.names = F, sep = '\t')

dtranscript = subset(d, V3 == 'transcript')
dtranscript_ann = sapply(as.vector(dtranscript$V9), function(x){
  # tmp = do.call(rbind, strsplit(strsplit(as.character(x), '\\;')[[1]], '\\='))
  tmp = do.call(rbind, strsplit(strsplit(as.character(x), '\\; ')[[1]], '\\ '))
  ttmp =as.vector(tmp[,2]); names(ttmp) = as.vector(tmp[,1])
  return(ttmp[c('gene_id', 'transcript_id', 'transcript_type', 'gene_name', 'transcript_name')])
})
dtranscript_ann = t(dtranscript_ann)
rownames(dtranscript_ann) = seq(1, nrow(dtranscript_ann))
dtranscript_ann = cbind(dtranscript[,c(1,4,5,7)], dtranscript_ann)
dtranscript_ann = dtranscript_ann[!duplicated(dtranscript_ann$transcript_id),]
dexon = subset(d, V3 == 'exon')
dexon_ann = sapply(as.vector(dexon$V9), function(x){
  tmp = do.call(rbind, strsplit(strsplit(gsub('  ', ' ', as.character(x)), '\\; ')[[1]], '\\ '))
  ttmp =as.vector(tmp[,2]); names(ttmp) = as.vector(tmp[,1])
  return(ttmp[c('transcript_id', 'exon_number')])
})
dexon_ann = t(dexon_ann)
rownames(dexon_ann) = seq(1, nrow(dexon_ann))
dexon_ann = cbind(dexon[,c(1, 4, 5, 7)], dexon_ann)

exonStart = tapply(dexon_ann[,2], factor(dexon_ann[,5]), function(x){
  paste0(x, collapse = ',')
})
exonEnd = tapply(dexon_ann[,3], factor(dexon_ann[,5]), function(x){
  paste0(x, collapse = ',')
})
exonNum = tapply(dexon_ann[,6], factor(dexon_ann[,5]), function(x){
  max(as.numeric(as.vector(x)))
})

gencode_lncRNA_ann = cbind(dtranscript_ann, exonStart[as.vector(dtranscript_ann$transcript_id)],
                           exonEnd[as.vector(dtranscript_ann$transcript_id)],
                           exonNum[as.vector(dtranscript_ann$transcript_id)])
gencode_lncRNA_ann$ID = do.call(rbind, strsplit(as.vector(gencode_lncRNA_ann$gene_id), '\\.'))[,1]
gencode_lncRNA_ann$length = abs(gencode_lncRNA_ann$V5-gencode_lncRNA_ann$V4)
colnames(gencode_lncRNA_ann) = c('chr', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'transcript_type', 'gene_name', 'transcript_name', 'exonStart', 'exonEnd', 'exonNum', 'ID', 'length')
# genes = do.call(rbind, strsplit(as.vector(gencode_lncRNA_ann$gene_id), '\\.'))[,1]
# genes2 = do.call(rbind, strsplit(as.vector(rownames(exp_mat)), '\\.'))[,1]
write.table(gencode_lncRNA_ann, 'gencode19release_lncRNA_feature.txt', quote = F, row.names = F, sep = '\t')
