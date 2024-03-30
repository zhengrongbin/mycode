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
    make_option(c("-c", "--count"), type = "character", default=TRUE,
              help="path for count matrix, should be csv file, unique gene name in the first columns"),
    make_option(c("-t", "--tpm"), type = "character", default=TRUE,
              help="path for TPM matrix, should be csv file, unique gene name in the first columns"),
    make_option(c('-m', "--meta_desgin"), type = 'character', default=TRUE,
              help='design matrix for specifying comparison, columns are samples, rows are comparisons with first column as comparison name, then 1 as treatment, 0 as control, keep empty for other samples'),
    make_option(c("-g", "--gene_ont"), type = "character", default=FALSE,
              help="True or False, True for doing gene ontology for up and down genes, False for no"),
    make_option(c("-e", "--gsea"), type = "character", default=FALSE,
              help="True or False, True for doing GSEA for the KEGG, False for no"),
    make_option(c("-s", "--species"), type = "character", default=FALSE,
              help="species given by hsa for human, and mmu for mouse, fish for zebrafish"),
    make_option(c("-p", "--prefix"), type = "character", default=FALSE,
              help="the prefix of output"),
    make_option(c("-o", "--protein_coding"), type = "character", default=FALSE,
              help="only protein coding, the gene annotation file, gene_name as column to notify the gene name"),
    make_option(c("-k", "--kegg"), type = "character", default=FALSE,
              help="the rds path for kegg gene list for fgsea enrichment"),
    make_option(c("-r", "--Replicate"), type = "character", default=FALSE,
              help="True or False, True for correcting replicates in DESeq2, make sure sample name ends with _rep1, _rep2, .etc")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

count_path <- opt$count
tpm_path <- opt$tpm
design_path <- opt$meta_desgin
protein_coding_path <- opt$protein_coding
prefix <- opt$prefix
species <- opt$species
go <- opt$gene_ont
gsea <- opt$gsea
kegg_path = opt$kegg
Replicate = opt$Replicate

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

    ## KEGG
    gene_map = bitr(geneset, fromType='SYMBOL', toType='ENTREZID', OrgDb=db, drop = TRUE)
    rownames(gene_map) = as.character(gene_map[,2])
    geneid = as.vector(gene_map[,2])
    print(head(geneid))
    print('KEGG')
    if (background == 'all'){
        kegg_res = enrichKEGG(gene = geneid, organism=species, keyType='kegg', pvalueCutoff = 1, qvalueCutoff = 1)
      }else{
        background_genes = as.vector(read.csv(background, header = F, sep = '\t')[,1])
        kegg_res = enrichKEGG(gene = geneid, organism=species, keyType='kegg', pvalueCutoff = 1, qvalueCutoff = 1, universe = background_genes)
     }
     GO[['KEGG']] = kegg_res
     return(list('gene_map'=gene_map, 'go'=GO))
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
    g = enrichplot::dotplot(GO[['KEGG']], showCategory = 20, title = paste0(label, '_KEGG'))
    print(g)
    dev.off()
}


## vocanoplot
vocano_plot <- function(mat, prefix){
    plot_df = subset(mat, (!is.na(padj)) & (!is.na(log2FoldChange)))
    # plot_df$diff = (abs(plot_df$log2FoldChange) > log2(1.5)) & (plot_df$padj < 0.01)
    # plot_df$diff = factor(plot_df$diff, levels = c('TRUE', 'FALSE'))
    plot_df$diff = 'Stable'
    plot_df$diff[(plot_df$log2FoldChange > log2(1.5)) & (plot_df$padj < 0.05)] = 'Up'
    plot_df$diff[(plot_df$log2FoldChange < -log2(1.5)) & (plot_df$padj < 0.05)] = 'Down'

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
cmat <- read.csv(count_path, row.names = 1)
if (tpm_path != NULL & tpm_path != FALSE){
tpm <- read.csv(tpm_path, row.names = 1)
}
design <- read.csv(design_path, row.names = 1)
## extract comparison
comparisons <- lapply(rownames(design), function(x){
    res = t(design)[,x]
    res[!is.na(res)]
})
names(comparisons) <- rownames(design)

### do PCA and hclust
print('++++ PCA and hcluster ++++')
if (tpm_path != NULL & tpm_path != FALSE){
    tpm.log = log10(tpm+1)
}else{
    tpm.log = log10(cmat+1)
}
## remove all zero
tpm.log = tpm.log[-which(apply(tpm.log, 1, function(x) all(x == 0))),]
var.3k = names(sort(apply(tpm.log, 1, var), decreasing = T)[1:3000])
tpm.log.scale.var = t(scale(t(tpm.log)))[var.3k,]

pc.cr=prcomp(t(tpm.log))
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

pdf(paste0(prefix, '_pca.pdf'), width = 5, height = 4.5)
print(g)
dev.off()

## hclust
pdf(paste0(prefix, '_hclust.pdf'), width = 5, height = 4.5)
plot(hclust(dist(t(tpm.log))))
dev.off()
## top variable gene heatmap
pdf(paste0(prefix, '_heatmap.pdf'), width = 5, height = 5.5)
Heatmap(tpm.log.scale.var, show_row_names = F, name='Scaled Exp', show_row_dend = T)
dev.off()

### do each DESeq2
diff_res = list(allgene=list(), protein_only=list())
for (comp in names(comparisons)) {
    if (file.exists(paste0(prefix, '_', comp, ".allgene.csv"))){
        diffexp_res <- read.csv(paste0(prefix, '_', comp, ".allgene.csv"), row.names = 1)
        diff_res[['allgene']][[comp]] <- diffexp_res
        vocano_plot(mat=diffexp_res, prefix = paste0(prefix, '_', comp, '_allgene'))
        next
    }
    if (Replicate == 'True'){
        samples <- names(comparisons[[comp]])
        cond = DataFrame('cond'=comparisons[[comp]], 'Replicate'=sapply(samples, function(x) tail(strsplit(x, '\\_')[[1]], 1)))
        tmp_count <- cmat[,samples]
        dds <- DESeqDataSetFromMatrix(tmp_count,
                                    cond, ~ Replicate+cond)
    }else{
        cond = DataFrame('cond'=comparisons[[comp]])
        samples <- names(comparisons[[comp]])
        tmp_count <- cmat[,samples]
        dds <- DESeqDataSetFromMatrix(tmp_count,
                                    cond, ~ cond)
    }
    dds <- DESeq(dds)
    ## save DEseq obj
    saveRDS(dds, file = paste0(prefix, '_', comp, '_deseq2_obj.rds'))
    diffexp_res <- as.data.frame(results(dds, contrast=list('cond')))
    ## plot
    vocano_plot(mat=diffexp_res, prefix = paste0(prefix, '_', comp, '_allgene'))
    ## write out
    write.csv(subset(diffexp_res, !is.na(log2FoldChange) & !is.na(padj)), file = paste0(prefix, '_', comp, ".allgene.csv"), quote = F)
    diff_res[['allgene']][[comp]] <- diffexp_res
    rm(diffexp_res)
}

##### only protein coding 
if (file.exists(protein_coding_path) == TRUE){
    print('++++ protein coding ++++++')
    gene_ann = read.csv(protein_coding_path)
    protein_coding = as.vector(gene_ann$gene_name)

    diff_res_protein = list()
    shared_protein_genes = intersect(protein_coding, rownames(cmat))
    for (comp in names(comparisons)) {
        if (file.exists(paste0(prefix, '_', comp, ".protein_gene.csv"))){
            pdiffexp_res <- read.csv(paste0(prefix, '_', comp, ".protein_gene.csv"), row.names = 1)
            diff_res[['protein_only']][[comp]] <- pdiffexp_res
            vocano_plot(mat=pdiffexp_res, prefix = paste0(prefix, '_', comp, '_protein_coding'))
            next
        }
        # cond = DataFrame('cond'=comparisons[[comp]])
        # samples <- names(comparisons[[comp]])
        # tmp_count <- cmat[shared_protein_genes,samples] ## only protein coding genes
        # dds <- DESeqDataSetFromMatrix(tmp_count,
        #                             cond, ~ cond)
        # dds <- DESeq(dds)
        # ## write out
        # pdiffexp_res = as.data.frame(results(dds))
        ddf = diff_res[['allgene']][[comp]]
        pdiffexp_res = subset(ddf, rownames(ddf) %in% shared_protein_genes)
        ## 
        vocano_plot(mat=pdiffexp_res, prefix = paste0(prefix, '_', comp, '_protein_coding'))

        write.csv(subset(pdiffexp_res, !is.na(log2FoldChange) & !is.na(padj)),, file = paste0(prefix, '_', comp, ".protein_gene.csv"), quote = F)
        diff_res[['protein_only']][[comp]] <- pdiffexp_res
        rm(pdiffexp_res)
    }
}


### functional analysis
## GO for up and down genes +++ only protein coding genes
if (go == 'True'){
    print('+++++ Gene Ontology ++++')
    for (comp in names(diff_res$protein_only)){
        if (length(diff_res$protein_only) > 0){
            diffexp = diff_res$protein_only[[comp]]
        }else{
            diffexp = diff_res$allgene[[comp]]
        }
        up = subset(diffexp, log2FoldChange >= log2(1.5) & padj <= 0.05)
        down = subset(diffexp, log2FoldChange <= -log2(1.5) & padj <= 0.05)
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
                if (k == 'KEGG'){
                    gene_map = up_go_res$gene_map
                    rownames(gene_map) = as.character(gene_map[,2])
                    res$geneID = sapply(as.vector(res$geneID), function(x){paste0(gene_map[strsplit(x, '\\/')[[1]], 'SYMBOL'], collapse = '/')})
                }
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
                if (k == 'KEGG'){
                    gene_map = down_go_res$gene_map
                    rownames(gene_map) = as.character(gene_map[,2])
                    res$geneID = sapply(as.vector(res$geneID), function(x){paste0(gene_map[strsplit(x, '\\/')[[1]], 'SYMBOL'], collapse = '/')})
                }
                write.csv(res, file = paste0(prefix, '_', comp, '.DOWN.', k, '.csv'), quote = F)
                rm(res)
            }
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
                       plot.title = element_text(hjust=0.5,vjust = 0.5,
                                                 margin = margin(l=100,r=50,t=10,b=10),face = "bold", colour = "black"))+
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
    saveRDS(fgseaRes, file = paste0(prefix, '_fgseaRes.rds'))
    return(plot_df)
}

### GSEA
if (gsea == "True"){
    print('+++++ GSEA')
    for (comp in names(diff_res$protein_only)){
        print(comp)
        diffexp <- diff_res$protein_only[[comp]]
        diffexp <- subset(diffexp, !is.na(log2FoldChange) & !is.na(padj))

        ### using log2FoldChange
        gageinput <- as.vector(diffexp$log2FoldChange)
        names(gageinput) <- as.vector(rownames(diffexp))
        gageinput <- gageinput[!grepl('^mt-|MT-', names(gageinput))]

        geneid <- bitr(names(gageinput), fromType='SYMBOL', toType='ENTREZID', OrgDb=db, drop = TRUE)
        rownames(geneid) <- as.character(geneid$ENTREZID)

        gageinput <- gageinput[as.vector(geneid$SYMBOL)]
        names(gageinput) <- as.vector(geneid$ENTREZID)

        gseainput = sort(gageinput, decreasing=TRUE)
        gseainput = gseainput[is.finite(gseainput)]
        print(species)
        print(length(gseainput))
        fullgsea <- gseKEGG(geneList = gseainput,
                        organism= species,
                        nPerm= 1000,
                        minGSSize= 50,
                        maxGSSize = 500,
                        pvalueCutoff = 0.99,
                        verbose= FALSE,
                        use_internal_data = FALSE)
        ## write out
        saveRDS(fullgsea, file = paste0(prefix, comp, '.gseaKEGG.rds'))
        pdf(paste0(prefix, comp, '.gseaKEGG.pdf'), width = 10, height = 8.5)
        g = dotplot(fullgsea, showCategory = 15, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
        print(g)
        g = ridgeplot(fullgsea) + labs(x = "Enrichment Distribution")
        print(g)
        dev.off()
        gsea_data = as.data.frame(fullgsea)
        gsea_data$core_enrichment = sapply(as.vector(gsea_data$core_enrichment), function(x){paste0(geneid[strsplit(x, '\\/')[[1]], 'SYMBOL'], collapse = '/')})
        gsea_data = gsea_data[order(-abs(gsea_data$NES)),]
        gsea_data[,2] = gsub(",", "", as.matrix(gsea_data[,2]))
        write.csv(gsea_data, file = paste0(prefix, '_', comp, '.gseaKEGG.csv'))
    }
    ## fgsea
    if (kegg_path != FALSE){
        kegg_gmt_list = readRDS(kegg_path) # GSEA reference

        for (comp in names(diff_res$protein_only)){
            print(comp)
            diffexp <- diff_res$protein_only[[comp]]
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
}

print('++++++ Finished +++++++')



