
library('optparse')

option_list <- list(
    make_option(c("-c", "--count"), type = "character", default=TRUE,
              help="path for count matrix, should be tab-separated file, unique peak cordinates in the first three columns"),
    make_option(c('-m', "--meta_desgin"), type = 'character', default=TRUE,
              help='design matrix for specifying comparison, columns are samples, rows are comparisons with first column as comparison name, then 1 as treatment, 0 as control, keep empty for other samples'),
    make_option(c('-s', "--stat_map"), type = 'character', default=FALSE,
              help='a table include information of mapping statistics from samtools of each bam'),
    make_option(c("-e", "--species"), type = "character", default=FALSE,
              help="species given by hsa for human, and mmu for mouse, fish for zebrafish"),
    make_option(c("-p", "--prefix"), type = "character", default=FALSE,
              help="the prefix of output"),
    make_option(c("-o", "--directory"), type = "character", default=FALSE,
              help="the directory of output files")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

count_path <- opt$count
design_path <- opt$meta_desgin
map_stat <- opt$stat_map
prefix <- opt$prefix
species <- opt$species
outdir <- opt$directory

# if (species == 'hsa'){
#     db = org.Hs.eg.db
# }else if (species == 'mmu'){
#     db = org.Mm.eg.db
#     species = 'mmu'
# }else{
#     db = org.Dr.eg.db
#     species = 'dre'
# }

if (prefix == FALSE){
    prefix = ''
}else{
    if (!grepl('_$', prefix)){
        prefix=paste0(prefix, '_')
    }
}

if (outdir == FALSE){
    outdir = './'
}


library('DESeq2')
library('ggplot2')
library('ggrepel')

## ======== read basic
print('++++ read in')
## ======== read count
cmat <- read.csv(count_path, sep = '\t')
rnames <- apply(cmat[,1:3], 1, function(x){gsub(' ', '', paste0(x, collapse = ':'))})
rownames(cmat) <- rnames
cmat <- cmat[,4:ncol(cmat)]

## total read count for each sample
if (map_stat == FALSE){
    sf = FALSE
}else{
    stats <- read.csv(map_stat, sep = '\t', row.names = 1) 
    ##
    sf <- stats[,1]
    names(sf) <- rownames(stats)
}

design <- read.csv(design_path, row.names = 1)
## extract comparison
comparisons <- lapply(rownames(design), function(x){
    res = t(design)[,x]
    res[!is.na(res)]
})
names(comparisons) <- rownames(design)

if (file.exists(paste0(outdir, '/', prefix, '.pca.pdf')) & file.exists(paste0(outdir, '/', prefix, '.hclust.pdf'))){
    print('')
}else{
    ### do PCA and hclust
    print('++++ PCA and hcluster ++++')
    cmat.log = log10(cmat+1)
    var.3k = names(sort(apply(cmat.log, 1, var), decreasing = T)[1:3000])
    cmat.log = cmat.log[var.3k,]

    pc.cr=prcomp(t(cmat.log))
    imp = summary(pc.cr)$importance

    ### plot pca
    plot_df = as.data.frame(pc.cr$x)
    plot_df$sname = rownames(plot_df)
    g = ggplot(data = plot_df, aes(x = PC1, y = PC2))+
            geom_point()+theme_bw()+
            theme(axis.text.x=element_text(size=10, colour = "black"),
            axis.text.y=element_text(size=10, colour = "black"),
            panel.border = element_blank(),axis.line = element_line(colour = "black"),
            text=element_text(size=14, colour = "black"),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust=0.5,vjust = 0.5,
                                    margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))+
            geom_text_repel(aes(label = sname), colour = 'black')+
            xlab(paste0('PC1 (', round(imp[2,'PC1']*100, 2), '%)'))+
            ylab(paste0('PC2 (', round(imp[2,'PC2']*100, 2), '%)'))

    pdf(paste0(outdir, '/', prefix, '.pca.pdf'), width = 4, height = 3.5)
    print(g)
    dev.off()

    ## hclust
    pdf(paste0(outdir, '/', prefix, '.hclust.pdf'), width = 4, height = 3.5)
    plot(hclust(dist(t(cmat.log))))
    dev.off()
}

## vocanoplot
vocano_plot <- function(mat, prefix){
    plot_df = subset(mat, (!is.na(padj)) & (!is.na(log2FoldChange)))
    # plot_df$diff = (abs(plot_df$log2FoldChange) > log2(1.5)) & (plot_df$padj < 0.01)
    # plot_df$diff = factor(plot_df$diff, levels = c('TRUE', 'FALSE'))
    plot_df$diff = 'Stable'
    plot_df$diff[(plot_df$log2FoldChange > log2(2)) & (plot_df$padj < 0.05)] = 'Up'
    plot_df$diff[(plot_df$log2FoldChange < -log2(2)) & (plot_df$padj < 0.05)] = 'Down'

    countDiff = table(plot_df$diff)

    plot_df$diff = factor(plot_df$diff, levels = c('Up', 'Stable', 'Down'))

    labels = rbind(head(plot_df[order(plot_df$stat),], 20), 
                    tail(plot_df[order(plot_df$stat),], 20))
    title = paste0('UP=', countDiff['Up'], '\nDOWN=', countDiff['Down'])

    labels$gene = rownames(labels)
    g = ggplot(data = plot_df,
            aes(x = log2FoldChange, y = -log10(padj), colour = diff))+
            geom_point()+theme_bw()+
            theme(axis.text.x=element_text(size=10, colour = "black"),
            axis.text.y=element_text(size=10, colour = "black"),
            panel.border = element_blank(),axis.line = element_line(colour = "black"),
            text=element_text(size=14, colour = "black"),
            legend.text=element_text(size=10),
            plot.title = element_text(colour = "black"))+
            scale_colour_manual(values=c('Up'='red', 'Stable'='grey', 'Down'='#1E90FF'))+
            labs(title = title)
            # geom_text_repel(data = labels, aes(x = log2FoldChange, y = -log10(padj), label = gene), colour = 'black')
    pdf(paste0(prefix, '_volcano.pdf'), width = 8)
    print(g)
    dev.off()
}

### do each DESeq2
diff_res = list()
for (comp in names(comparisons)) {
    if (!file.exists(paste0(outdir, '/', prefix, comp, '_deseq2_obj.rds'))){
        cond = DataFrame('cond'=comparisons[[comp]])
        samples <- names(comparisons[[comp]])
        tmp_count <- cmat[,samples]
        dds <- DESeqDataSetFromMatrix(tmp_count,
                                    cond, ~ cond)
        dds <- dds[rowSums(counts(dds))>0,]
        if (sf == FALSE){
            dds <- DESeq(dds)
        }else{
            sizeFactors(dds) <- sf[samples] / mean(sf[samples])
            dds <- estimateDispersions(dds)
            dds <- nbinomWaldTest(dds)
        }  
        ## save DEseq obj
        saveRDS(dds, file = paste0(outdir, '/', prefix, comp, '_deseq2_obj.rds'))
    }else{
        dds <- readRDS(paste0(outdir, '/', prefix, comp, '_deseq2_obj.rds'))
        pdf(paste0(outdir,'/',prefix, comp, '_MA.pdf'))
        plotMA(results(dds))
        dev.off()
        diffpeak_res <- as.data.frame(results(dds))
        ## plot
        vocano_plot(mat=diffpeak_res, prefix = paste0(outdir,'/', prefix, comp))
        ## write out
        write.csv(subset(diffpeak_res, !is.na(log2FoldChange) & !is.na(padj)), file = paste0(outdir, '/', prefix, comp, '_diffres.csv'), quote = F)
        diff_res[[comp]] <- subset(diffpeak_res, !is.na(log2FoldChange) & !is.na(padj))
        rm(diffpeak_res)
    }
}

## further analysis of diff peaks
for (comp in names(diff_res)){
    print(paste0('+++ diff peaks for', comp))
    up = subset(diff_res[[comp]], log2FoldChange > 1 & padj < 0.05)
    up_peaks = do.call(rbind, strsplit(rownames(up), '\\:'))
    write.table(up_peaks, paste0(outdir,'/', prefix, comp, '_UPpeak.bed'), sep = '\t', quote = F, row.names = F, col.names = F)

    down = subset(diff_res[[comp]], log2FoldChange < -1 & padj < 0.05)
    down_peaks = do.call(rbind, strsplit(rownames(down), '\\:'))
    write.table(down_peaks, paste0(outdir,'/',prefix, comp, '_DOWNpeak.bed'), sep = '\t', quote = F, row.names = F, col.names = F)

    shared = subset(diff_res[[comp]], abs(log2FoldChange) < 1 | padj > 0.05)
    shared_peaks = do.call(rbind, strsplit(rownames(shared), '\\:'))
    write.table(shared_peaks, paste0(outdir,'/',prefix, comp, '_SHAREDpeak.bed'), sep = '\t', quote = F, row.names = F, col.names = F)
}





