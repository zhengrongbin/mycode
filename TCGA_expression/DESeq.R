args = commandArgs(T)
library('DESeq2')
library('ggplot2')
library('ggrepel')

count_path = args[1] # count matrix
design_path = args[2] # design matrix
prefix = args[3] # prefix of outputs

static = './static/'
protein_coding_ann_path = paste0(static, 'gene_annotation_GENCODE.v22.csv')

config = list(
    "count_path" = count_path,
    "design_path" = design_path,
    'prefix' = prefix,
    'protein_coding' = protein_coding_ann_path
)
print(config)

# mesg
msg = function(cnt){
    system(paste0('echo "++++++ ', cnt, '"'))
}

# read in count file
read_in = function(count_path, design_path){
    msg(count_path)
    msg(design_path)
    counts = read.csv(count_path, row.names = 1)
    design = read.csv(design_path, row.names = 1)
    if (ncol(counts) < 4){
        msg('too few samples for DESeq2, less than 4 was found.')
        quit()
    }
    ## whether all the count column names are in design table
    design_match <- length(intersect(as.vector(rownames(design)), 
                                as.vector(colnames(counts)))) == ncol(counts)
    if (design_match != T){
        msg('design table does match to count columns')
    }
    design_match <- intersect(as.vector(rownames(design)), 
                                as.vector(colnames(counts)))
    design = design[design_match,,drop=FALSE]
    counts = counts[,design_match,drop=FALSE]
    return(list('count'=counts, 'desgin'=design))
}


deseq = function(counts, design, prefix){
    msg('running DESeq2')
    # run DESeq2 and get result of each comparison
    cond <- DataFrame("cond"=gsub(' ', '_', as.vector(design[,1])))
    dds <- DESeqDataSetFromMatrix(counts,
                                cond, ~ cond)
    dds <- DESeq(dds)
    diff_res <- list()
    for (compare in resultsNames(dds)){
        if (compare != 'Intercept'){
            res <- results(dds, contrast = list(compare))
            msg(paste0('output diff exp of ', compare))
            write.csv(res, file = paste0(prefix, '_', compare, '.csv'), quote = F)
            diff_res[[compare]] <- as.data.frame(res)
        }
    }
    return(diff_res)
}

## vocanoplot
vocano_plot <- function(mat, prefix){
    plot_df = subset(mat, (!is.na(padj)) & (!is.na(log2FoldChange)))
    plot_df$diff = (abs(plot_df$log2FoldChange) > log2(1.5)) & (plot_df$padj < 0.01)
    plot_df$diff = factor(plot_df$diff, levels = c('TRUE', 'FALSE'))
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
            scale_colour_manual(values=c('TRUE'='red', 'FALSE'='grey'))+
            geom_text_repel(data = labels, aes(x = log2FoldChange, y = -log10(padj), label = gene), colour = 'black')
    pdf(paste0(prefix, '_volcano.pdf'), width = 8)
    print(g)
    dev.off()
}

run = function(config){
    msg(config['protein_coding'])
    protein_coding <- read.csv(config[['protein_coding']])
    print('===')
    tmp <- read_in(config[['count_path']], config[['design_path']])
    diff_exp <- deseq(counts = tmp[['count']], design = tmp[['desgin']], prefix=config[['prefix']])
    ## plot diff exp in volcano
    for (compare in names(diff_exp)){
        msg('plotting volcano')
        vocano_plot(subset(diff_exp[[compare]], rownames(diff_exp[[compare]]) %in% protein_coding$gene_name),
                     paste0(prefix, '_coding_gene_', compare)) # only coding gene
        vocano_plot(subset(diff_exp[[compare]], !rownames(diff_exp[[compare]]) %in% protein_coding$gene_name),
                     paste0(prefix, '_non-coding_gene_', compare)) # only coding gene
    }
}

## excute
run(config)
