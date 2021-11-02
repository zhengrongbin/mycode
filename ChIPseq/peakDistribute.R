library(ChIPseeker)
library(clusterProfiler)

dirpath = './'
species = 'mm10'

if (species == 'mm10'){
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
}else if (species == 'hg38'){
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
}else if (species == 'hg19'){
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
}else{
	library(TxDb.Mmusculus.UCSC.mm9.knownGene)
	txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
}

setwd(dirpath)
files = list.files(pattern = 'summits.bed$', path = './', full.name = T)
names(files) = gsub('.rep1_sorted_summits.bed', '', basename(files))#sapply(files, function(x) strsplit(basename(x), '\\.')[[1]][1])

# peaj_distribution_summary = list()
# pdf('peak_pie.pdf')
# files = list.files(path = './', pattern = '.bed$')
# for (file in files){
#     peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000),
#                          TxDb=txdb)
#     cond = strsplit(file, '\\.')[[1]][1]
#     g = plotAnnoPie(peakAnno)
#     text(0.1,g$rect$top,cond)
#     peaj_distribution_summary[[cond]] = peakAnno 
# }
# dev.off()

if (length(files) > 0){
    peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                           tssRegion=c(-3000, 3000), verbose=FALSE)
    print(peakAnnoList)
    # names(peakAnnoList) = gsub('.rep1_peaks.narrowPeak|.//', '', files)
    print('+++++ narrowPeak')
    pdf(paste0(dirpath, '_peakSummit_distribute.pdf'))
    p = plotAnnoBar(peakAnnoList)
    print(p)
    dev.off()
    saveRDS(peakAnnoList, file = paste0(dirpath, 'peakSummit_distribute.rds'))
    write.csv(do.call(rbind, lapply(peakAnnoList, function(x){x@annoStat})), file = 'peakSummit_distribute.csv', quote = F)
}

## ========= up peak
files = list.files(pattern = '.up.bed$', path = './', full.name = T)
names(files) = sapply(files, function(x) strsplit(basename(x), '\\.')[[1]][1])

if (length(files) > 0){
    peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                           tssRegion=c(-3000, 3000), verbose=FALSE)
    print('+++++ up peak')
    pdf(paste0(dirpath, '_diff_up_peak_distribute.pdf'))
    p1 = plotAnnoBar(peakAnnoList)
    print(p1)
    dev.off()
    saveRDS(peakAnnoList, file = paste0(dirpath, '_diff_up_peak_distribute.rds'))
}

## ========= down peak
files = list.files(pattern = '.down.bed$', path = './', full.name = T)
names(files) = sapply(files, function(x) strsplit(basename(x), '\\.')[[1]][1])

if (length(files) > 0){
    peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                           tssRegion=c(-3000, 3000), verbose=FALSE)
    print('+++++ down peak')
    pdf(paste0(dirpath, '_diff_down_peak_distribute.pdf'))
    p2 = plotAnnoBar(peakAnnoList)
    print(p2)
    dev.off()
    saveRDS(peakAnnoList, file = paste0(dirpath, '_diff_down_peak_distribute.rds'))
}
print(paste0('+++ finished: ', dirpath))