## a DESeq2 differential expression analysis workfow
## find differentially expressed genes and followed with functional enrichment

library('DESeq2')
library('fgsea')
library('ggplot2')
library('ggrepel')
library('optparse')
library('clusterProfiler')
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library('org.Dr.eg.db')
library('ComplexHeatmap')

option_list <- list(
    make_option(c("-i", "--input_diff"), type = "character", default=FALSE,
              help="File path of DESeq2 differential expression table"),
    make_option(c("-f", "--fc_cut"), default=1.5,
              help="Fold Change Cutoff"),
    make_option(c("-q", "--padj_cut"), default=0.05,
              help="adjusted p value Cutoff"),
    make_option(c("-g", "--gene_ont"), type = "character", default=FALSE,
              help="True or False, True for doing gene ontology for up and down genes, False for no"),
    make_option(c("-e", "--gsea"), type = "character", default=FALSE,
              help="True or False, True for doing GSEA for the KEGG, False for no"),
    make_option(c("-s", "--species"), type = "character", default=FALSE,
              help="species given by hsa for human, and mmu for mouse, fish for zebrafish"),
    make_option(c("-p", "--prefix"), type = "character", default=FALSE,
              help="the prefix of output"),
    make_option(c("-k", "--kegg"), type = "character", default=FALSE,
              help="the rds path for kegg gene list for fgsea enrichment"),
    make_option(c("-l", "--label"), type = "character", default=FALSE,
              help="label for the comparison")
    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input_diff
prefix <- opt$prefix
species <- opt$species
go <- opt$gene_ont
gsea <- opt$gsea
kegg_path = opt$kegg
fc_cut <- as.numeric(opt$fc_cut)
padj_cut = as.numeric(opt$padj_cut)
comp = opt$label


if (species == 'hsa'){
    db = org.Hs.eg.db
}else if (species == 'mmu'){
    db = org.Mm.eg.db
    species = 'mmu'
}else{
    db = org.Dr.eg.db
    species = 'dre'
}

if (prefix == FALSE){
    prefix = './'
}

## ======== functions =======
enrichment_analysis <- function(geneset, background = 'all', species = 'hsa'){
    GO = list()
    for (k in c('BP', 'MF', 'CC')){
      print(k)
      if (background == 'all'){
        GO[[k]] = enrichGO(gene = geneset, OrgDb = db, pvalueCutoff = 1, qvalueCutoff = 1, ont = k, keyType = 'SYMBOL')
      }else{
        background_genes = as.vector(read.csv(background, header = F, sep = '\t')[,1])
        GO[[k]] = enrichGO(gene = geneset, OrgDb = db, pvalueCutoff = 1, qvalueCutoff = 1, ont = k, universe = background_genes, keyType = 'SYMBOL')
      }
    }
     return(list('go'=GO))
}

## plot GO
plot_go <- function(GO, label, prefix){
    pdf(paste0(prefix, label, '.pdf'), width = 10, height = 6.5)
    g = enrichplot::dotplot(GO[['MF']], showCategory = 20, title = paste0(label, '_MF'))
    print(g)
    g = enrichplot::dotplot(GO[['BP']], showCategory = 20, title = paste0(label, '_BP'))
    print(g)
    g = enrichplot::dotplot(GO[['CC']], showCategory = 20, title = paste0(label, '_CC'))
    print(g)
    # g = enrichplot::dotplot(GO[['KEGG']], showCategory = 20, title = paste0(label, '_KEGG'))
    # print(g)
    dev.off()
}


## vocanoplot
vocano_plot <- function(mat, prefix, fc_cut = 1.5, padj_cut = 0.05){
    plot_df = subset(mat, (!is.na(padj)) & (!is.na(log2FoldChange)))
    # plot_df$diff = (abs(plot_df$log2FoldChange) > log2(1.5)) & (plot_df$padj < 0.01)
    # plot_df$diff = factor(plot_df$diff, levels = c('TRUE', 'FALSE'))
    plot_df$diff = 'Stable'
    plot_df$diff[(plot_df$log2FoldChange > log2(fc_cut)) & (plot_df$padj < padj_cut)] = 'Up'
    plot_df$diff[(plot_df$log2FoldChange < -log2(fc_cut)) & (plot_df$padj < padj_cut)] = 'Down'

    plot_df$diff = factor(plot_df$diff, levels = c('Up', 'Stable', 'Down'))

    labels = rbind(head(plot_df[order(plot_df$stat),], 20), 
                    tail(plot_df[order(plot_df$stat),], 20))

    labels$gene = rownames(labels)
    g = ggplot(data = plot_df,
            aes(x = log2FoldChange, y = -log10(padj), colour = diff))+
            geom_point()+theme_bw()+
            theme(axis.text.x=element_text(size=10, colour = "black"),
            axis.text.y=element_text(size=10, colour = "black"),
            panel.border = element_blank(),axis.line = element_line(colour = "black"),
            text=element_text(size=14, colour = "black"),
            legend.text=element_text(size=10),
            plot.title = element_text(hjust=0.5,vjust = 0.5,
                                    margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))+
            scale_colour_manual(values=c('Up'='red', 'Stable'='grey', 'Down'='#1E90FF'))+
            labs(title=paste0('Up=', nrow(plot_df[plot_df[,'diff'] == 'Up',]), '\nDown=',
                             nrow(plot_df[plot_df[,'diff'] == 'Down',])))+
            geom_text_repel(data = labels, aes(x = log2FoldChange, y = -log10(padj), label = gene), colour = 'black')
    pdf(paste0(prefix, '_volcano.pdf'), width = 8)
    print(g)
    dev.off()
}
## ======== read basic
print('++++ read in')

## ======== read count
diffexp <- read.csv(input_file, row.names = 1)

## volcano plot
vocano_plot(mat=diffexp, prefix = paste0(prefix, '_', comp), fc_cut = fc_cut, padj_cut = padj_cut)


### functional analysis
## GO for up and down genes +++ only protein coding genes
if (go == 'True'){
    print('+++++ Gene Ontology ++++')
    up = subset(diffexp, log2FoldChange >= log2(fc_cut) & padj < padj_cut)
    down = subset(diffexp, log2FoldChange <= -log2(fc_cut) & padj < padj_cut)
    if (nrow(up) != 0){
        if (file.exists(file = paste0(prefix, '_', comp, ".GO_UP.rds"))){
            up_go_res = readRDS(file = paste0(prefix, '_', comp, ".GO_UP.rds"))
            plot_go(GO=up_go_res, label=paste0(comp, '_up_GO'), prefix = prefix)
            next
        }
        ### write out
        write.csv(up, file = paste0(prefix, '_', comp, ".protein_gene.UP.csv"), quote = F)
        ### functional enrichment
        # ===== up ======
        up_go_res = enrichment_analysis(geneset = rownames(up)[!grepl('^mt-|MT-', rownames(up))], background='all', species = species)
        ## plot
        saveRDS(up_go_res$go, file = paste0(prefix, '_', comp, ".GO_UP.rds"))
        plot_go(GO=up_go_res$go, label=paste0(comp, '_up_GO'), prefix = prefix)
        ### write out
        ## KEGG entrezid to gene name
        for (k in names(up_go_res[['go']])){
            res = as.data.frame(up_go_res[['go']][[k]])
            # if (k == 'KEGG'){
            #     gene_map = up_go_res$gene_map
            #     rownames(gene_map) = as.character(gene_map[,2])
            #     res$geneID = sapply(as.vector(res$geneID), function(x){paste0(gene_map[strsplit(x, '\\/')[[1]], 'SYMBOL'], collapse = '/')})
            # }
            write.csv(res, file = paste0(prefix, '_', comp, '.UP.', k, '.csv'), quote = F)
            rm(res)
        }
    }

    # ===== down ======
    if (nrow(down) != 0){
        if (file.exists(file = paste0(prefix, '_', comp, ".GO_DOWN.rds"))){
            up_go_res = readRDS(file = paste0(prefix, '_', comp, ".GO_DOWN.rds"))
            plot_go(GO=up_go_res, label=paste0(comp, '_down_GO'), prefix = prefix)
            next
        }
        write.csv(down, file = paste0(prefix, '_', comp, ".protein_gene.DOWN.csv"), quote = F)

        down_go_res = enrichment_analysis(geneset = rownames(down)[!grepl('^mt-|MT-', rownames(down))], background='all', species = species)
        saveRDS(down_go_res$go, file = paste0(prefix, '_', comp, ".GO_DOWN.rds"))
        plot_go(GO=down_go_res$go, label=paste0(comp, '_down_GO'), prefix = prefix)
        ### write out
        ## KEGG entrezid to gene name
        for (k in names(down_go_res[['go']])){
            res = as.data.frame(down_go_res[['go']][[k]])
            # if (k == 'KEGG'){
            #     gene_map = down_go_res$gene_map
            #     rownames(gene_map) = as.character(gene_map[,2])
            #     res$geneID = sapply(as.vector(res$geneID), function(x){paste0(gene_map[strsplit(x, '\\/')[[1]], 'SYMBOL'], collapse = '/')})
            # }
            write.csv(res, file = paste0(prefix, '_', comp, '.DOWN.', k, '.csv'), quote = F)
            rm(res)
        }
    }
    
}

fgsea_fuc <- function(rank_genes, kegg_gmt_list, prefix){
    set.seed(1234)
    fgseaRes <- fgsea(kegg_gmt_list, rank_genes[!is.na(rank_genes)],
                         minSize  = 15,
                         maxSize  = 500)

    fgseaRes$pathway = gsub('KEGG_', '', fgseaRes$pathway)
    fgseaRes = fgseaRes[order(-fgseaRes$NES),]
    plot_df = subset(fgseaRes, padj <= 0.1)
    saveRDS(fgseaRes, file = paste0(prefix, '_fgseaRes.rds'))
    g = ggplot(plot_df, aes(x = NES, y = reorder(pathway, NES), colour=padj, size = size))+
        geom_count()+
        geom_point(shape = 1, color = 'black')+
        xlab('Normalized Enrichment Score')+ylab('')+labs(title = 'Gene Set Enrichment Analysis')+
        theme_bw()+theme(axis.text.x=element_text(size=10, face = 'bold', colour = "black",
                                                angle = 45, hjust = 1, vjust = 1),
                       axis.text.y=element_text(size=8, colour = "black"),
                       panel.border = element_blank(),axis.line = element_line(colour = "black"),
                       text=element_text(size=10, colour = "black"),
                       legend.text=element_text(size=10),
                       plot.title = element_text(hjust=0.5,vjust = 0.5,margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))+
        scale_colour_gradient2(high = 'grey', low = 'red', limits = c(NA, 0.2), midpoint = 0.1)+
        geom_vline(xintercept = 0, color = 'grey')
    # save out gsea dot plot
    pdf(paste0(prefix, '_fgsea.pdf'), width = 10, height = 4 + (nrow(plot_df) * 0.1))
    print(g)
    dev.off()
    ## write out gsea result table
    res = as.matrix(fgseaRes)
    res[,'leadingEdge'] = gsub(',', ';', as.vector(res[,'leadingEdge']))
    write.table(res, file = paste0(prefix, '_fgseaRes.tsv'), quote = F, row.names = F, sep = '\t')   
    return(plot_df)
}

### GSEA
if (gsea == "True"){
    print('+++++ fGSEA')
    ## fgsea
    if (kegg_path != FALSE){
        kegg_gmt_list = readRDS(kegg_path) # GSEA reference

        diffexp <- subset(diffexp, !is.na(log2FoldChange) & !is.na(padj))

        ### using log2FoldChange
        gageinput <- as.vector(diffexp$log2FoldChange)
        names(gageinput) <- as.vector(rownames(diffexp))
        gageinput <- gageinput[!grepl('^mt-|MT-', names(gageinput))]

        fgsea_res = fgsea_fuc(rank_genes=gageinput,
                        kegg_gmt_list=kegg_gmt_list,
                        prefix=paste0(prefix, '_', comp))
    }
}

print('++++++ Finished +++++++')



