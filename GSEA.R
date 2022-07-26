
#!/usr/bin/env R
#description     :GSEA analysis
#author          :Rongbin Zheng
#date            :2020-08-29
#usage           :python prepare_data.py


args = commandArgs(T)
library('fgsea')
library('ggplot2')
library('ggrepel')
library('optparse')
# inputpath = args[1] #  input file where first column is gene name, other columns are ranking scores
# prefix = args[2] # prefix of outputs
# if (length(args) == 3){
#     column_sp = args[3]
# }else{
#     column_sp = NA
# }

# if (length(args) == 4){
#   kegg_path = args[4]
# }else{
#   static = './static/'
#   kegg_path = paste0(static, 'c2.cp.kegg.v7.1.symbols.gmt.rds')
# }

# if (length(args) == 0){
#     print('Rscript GSEA.R inputpath prefix [column_specified]')
#     quit()
# }
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default=FALSE,
              help="input file where first column is gene name, other columns are ranking scores"),
    make_option(c("-p", "--prefix"), type = "character", default=FALSE,
              help="the prefix of output"),
    make_option(c("-c", "--column"), type = "character", default=FALSE,
              help="specify a column in -i for ranking for GSEA"),
    make_option(c('-k', '--kegg'), type = 'character', default=FALSE,
              help = 'specify a KEGG reference, rds file with a list')
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

inputpath = opt$input
prefix = opt$prefix
column_sp = opt$column
kegg_path = opt$kegg
if (kegg_path == FALSE){
  static = '/Users/rongbinzheng/Documents/Workplace/Code/TCGA_Expression/static/'
  kegg_path = paste0(static, 'c2.cp.kegg.v7.1.symbols.gmt.rds')
  print(kegg_path)
}

config = list(
    "input_path" = inputpath,
    'prefix' = prefix,
    'kegg_path' = kegg_path,
    'column_sp' = column_sp
)
print(config)
# mesg
msg = function(cnt){
    system(paste0('echo "++++++ ', cnt, '"'))
}


gseaScores = function (geneList, geneSet, exponent = 1, fortify = FALSE) 
         {
         geneSet <- intersect(geneSet, names(geneList))
         N <- length(geneList)
         Nh <- length(geneSet)
         Phit <- Pmiss <- numeric(N)
         hits <- names(geneList) %in% geneSet
         Phit[hits] <- abs(geneList[hits])^exponent
         NR <- sum(Phit)
         Phit <- cumsum(Phit/NR)
         Pmiss[!hits] <- 1/(N - Nh)
         Pmiss <- cumsum(Pmiss)
         runningES <- Phit - Pmiss
         max.ES <- max(runningES)
         min.ES <- min(runningES)
         if (abs(max.ES) > abs(min.ES)) {
         ES <- max.ES
         }
         else {
         ES <- min.ES
         }
         df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
         position = as.integer(hits))
         if (fortify == TRUE) {
         return(df)
         }
         df$gene = names(geneList)
         res <- list(ES = ES, runningES = df)
         return(res)
}


gsInfo.new = function(geneList, geneSet, geneSetID, exponent=1){
         df.c = NULL
         # print(head(geneList))
         for (geneSetId in geneSetID){
            # print(names(geneSet))
           df <- gseaScores(geneList, geneSet[[geneSetId]], exponent, fortify=TRUE)
           df$ymin <- 0
           df$ymax <- 0
           pos <- df$position == 1
           h <- diff(range(df$runningScore))/20
           df$ymin[pos] <- -h
           df$ymax[pos] <- h
           df$geneList <- geneList

           df$Description <- geneSetId
           df.c = rbind(df.c, df)
         }
         return(df.c)
}
enrich_plot_custom = function(geneList, geneSet, geneSetID, pre_ranked = FALSE,
    title = "", color = "green", base_size = 11,
    rel_heights = c(1.5, 0.5, 1), subplots = 1:3, ylimit = c(),
    pvalue_table = FALSE, ES_geom = "line", col_low='Blues', col_high='Reds',
    alpha = .2)
    ## geneList -> gene score, such as foldchange, gene named
    ## geneSet -> a list of term - gene 
    ## geneSetID -> the term want to plot
    {
    if ( pre_ranked == FALSE ){
      geneList = sort(geneList,decreasing = T)
    }
    gsdata <- gsInfo.new(geneList = geneList, geneSet = geneSet,
              geneSetID = geneSetID)

    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(colour = "grey92"), 
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    #panel.grid.major = element_blank(), # get rid of major grid
    #panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")# get rid of legend panel bg 
          ) +
    scale_x_continuous(expand = c(0, 0))
    ## line
    p <- p + geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey')
    
    if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
    size = 1)
    } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
    size = 1, data = subset(gsdata, position == 1))
    }
    p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
    legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
    p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
    plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
    unit = "cm"))
    i <- 0
    for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
    term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
    ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) + theme(legend.position = "none", 
    plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
    axis.text = element_blank(), axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
    0))
    if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
    inv <- inv + 1
    col <- c(rev(RColorBrewer::brewer.pal(5, col_low)), RColorBrewer::brewer.pal(5, col_high))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
    ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
    alpha = 0.9, inherit.aes = FALSE)
    }
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
    y = ~y, yend = 0), color = "grey")
    p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
    l = 0.2, unit = "cm"))
    if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
        if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        }else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }
    ## add area color
    ## data frame of plotting from ggplot
    df = ggplot_build(p.res)$data[[2]]
    ## color set
    color_set = unique(df[,'colour'])
    names(color_set) = unique(gsdata[,'Description'])

    for (k in names(color_set)){
        p.res <- p.res + geom_area(data = subset(gsdata, Description == k), 
                      aes(x=x, y=runningScore), fill = color_set[k], alpha = .2)

    }
    ##
    if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
    xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
    0.95), ymin = quantile(p.res$data$runningScore, 
    0.75), ymax = quantile(p.res$data$runningScore, 
    0.9))
    }
    ## fgsea
    fgres = as.data.frame(fgsea(geneSet[geneSetID], geneList))
    rownames(fgres) = as.vector(fgres[,'pathway'])
    NES = fgres[geneSetID, 'NES']
    padj = fgres[geneSetID, 'padj']
    ## label
    p.res = p.res+labs(title = paste0(title, '\n', 'NES = ', round(NES, 2), ', padj = ', format(padj, digits = 2, scientific = T)))
    
    plotlist <- list(p.res, p2)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
    axis.ticks.x = element_line(), axis.text.x = element_text())
    if (length(subplots) == 1) {
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
    r = 0.2, b = 0.2, l = 0.2, unit = "cm")))}
    if (length(rel_heights) > length(subplots)) {
    rel_heights <- rel_heights[subplots]}
    g = cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
    return(g)
}

plotEnrichment <- function(pathway, pname, stats, gseaParam = 1, ticksSize = 0.2){
  stats <- stats[!is.na(stats)]
  breaks <- length(stats[stats>0])
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- statsAdj[!is.na(statsAdj)]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  gstyle = theme_bw()+theme(axis.text.x=element_text(size=8, colour = "black"),
                            text=element_text(size=7, colour = "black"),
                            legend.title = element_blank(),legend.text=element_text(size=10),
                            plot.title = element_text(hjust=0.5,vjust = 0.5, 
                                                      margin = margin(l=100,r=50,t=10,b=10),
                                                      face = "bold", colour = "black"))
  midstyle = theme(panel.grid = element_blank(),
                   axis.ticks.x=element_blank(), 
                   axis.text.x=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.y=element_blank(), plot.margin = unit(c(-0.7,1,-0.7,1), 'lines'))
  
  g1 = ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = "green", size = .5) +  
    geom_line(color = 'green', size = 0.7)+gstyle+
    geom_hline(yintercept = 0)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(), plot.margin = unit(c(1,1,-0.7,1), 'lines'))+xlab('')+
    labs(title = paste0('Enrichment Plot: \n', pname))+ylab('Enrichment Score (ES)')
  g2 = ggplot()+geom_segment(data = data.frame(x = pathway), 
                             mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) +
    gstyle+midstyle+ylab('')+xlab('')
    
  rank_mat = data.frame('ranks' = 0:(length(ord)-1), 'score'=stats[ord])
  rank_mat$h = 1
  g3 = ggplot(data = rank_mat, aes(x = ranks, y = h, fill = ranks))+geom_tile()+gstyle+midstyle+
    theme(legend.position='none')+
    xlab('')+ylab('')+
    scale_fill_gradient2(high = 'blue', low = 'red', mid = 'white', midpoint = breaks)
  
  g4 = ggplot(data = rank_mat[seq(1, nrow(rank_mat), 100),], aes(x = ranks, y = score))+
    geom_line(color = 'grey')+geom_area(fill = "lightgrey")+
    gstyle+theme(plot.margin = unit(c(-0.7,0,2,0), "lines"))+
    geom_vline(xintercept = breaks, colour = 'grey', linetype = 'dashed')+
    geom_text(x = breaks, y = 0, label = paste0('Zero cross at ', breaks), size = 2)+
    xlab('Rank in Ordered Dataset')+ylab('Ranked list metric (Preranked)')+
    geom_text(x = breaks + 50, y = min(rank_mat$score), label = 'negatively correlated', color = 'blue',
              hjust = 0, size = 2)+
    geom_text(x = 50, y = max(rank_mat$score), label = 'positively correlated', color = 'red',
              hjust = 0, size = 2)
  
  g = cowplot::plot_grid(g1, g2, g3, g4, align = "v", ncol = 1, rel_heights = c(0.5, 0.12, 0.03, 0.35))
  
  return(g)
}
                   


gsea <- function(rank_genes, kegg_gmt_list, prefix){
    set.seed(1234)
    fgseaRes <- fgsea(pathways = kegg_gmt_list, 
                         stats    = rank_genes[!is.na(rank_genes)],
                         minSize  = 15,
                         maxSize  = 500, nperm = 1000)

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
    msg('draw GSEA plot')
    pdf(paste0(prefix, '_GSEA.pdf'), width = 10, height = 4 + (nrow(plot_df) * 0.1))
    print(g)
    dev.off()
    ## write out gsea result table
    res = as.matrix(fgseaRes)
    res[,'leadingEdge'] = gsub(',', ';', as.vector(res[,'leadingEdge']))
    msg('save GSEA result')
    write.table(res, file = paste0(prefix, '_fgseaRes.tsv'), quote = F, row.names = F, sep = '\t')   
    saveRDS(fgseaRes, file = paste0(prefix, '_fgseaRes.rds'))
    return(plot_df)
}

plot_gsea_enrichment <- function(rank_genes, pathways, prefix, kegg_gmt_list){
  pdf(paste0(prefix,'_good_GSEA.pdf'), width = 4, height = 4.5)
  for (pname in pathways){
    # pname = paste0('KEGG_', pname)
    msg(paste0('GSEA ploting: ', pname))
    #pdf(paste0(prefix, gsub('\\/', '', pname), '_GSEA.pdf'), width = 4, height = 4.5)
    # g = plotEnrichment(kegg_gmt_list[[pname]], pname,
    #                    rank_genes, gseaParam = 1, ticksSize = 0.2)
    # print(head(rank_genes))
    g = enrich_plot_custom(geneList=rank_genes, geneSet=kegg_gmt_list, geneSetID=pname, pre_ranked = FALSE,
                                  title = pname, color = "green", base_size = 11,
                                  rel_heights = c(1.5, 0.5, 1), subplots = 1:3, ylimit = c(),
                                  pvalue_table = FALSE, ES_geom = "line", col_low='Blues', col_high='Reds',
                                  alpha = .2)
    print(g)
    #dev.off()
  }
  dev.off()
}

run = function(config){
    kegg_gmt_list = readRDS(config[["kegg_path"]]) # GSEA reference
    msg('read in diff exp matrix')
    # diff_exp = read.csv(deseq_path, row.names = 1)
    # rank_genes = diff_exp$stat # gsea using DESeq statistics
    inputFile = read.csv(inputpath, row.names = 1, check.names = F)
    if (config[['column_sp']] == FALSE){
        msg('no specified column, run all')
        idcs = colnames(inputFile)
        print(idcs)
        for (idc in idcs){
            rank_genes <- as.vector(inputFile[,idc])
            names(rank_genes) <- rownames(inputFile)
            msg(paste0('++++++ running GSEA for ', idc))
            good_pathway = gsea(rank_genes, kegg_gmt_list, paste0(config[['prefix']], '_', gsub('/', '', idc)))
            # plot_gsea_enrichment(rank_genes, as.vector(good_pathway$pathway), 
            #                      paste0(config[['prefix']], '_', idc), kegg_gmt_list)
        }
    } else {
        msg('specified a column')
        idc = config[['column_sp']]
        rank_genes <- as.vector(inputFile[,idc])
        names(rank_genes) <- rownames(inputFile)
        rank_genes <- rank_genes[!is.na(rank_genes)]
        msg(paste0('++++++ running GSEA for ', idc))
        good_pathway = gsea(rank_genes, kegg_gmt_list, paste0(config[['prefix']], '_', idc))
        plot_gsea_enrichment(rank_genes, as.vector(good_pathway$pathway),
                              paste0(config[['prefix']], '_', idc), kegg_gmt_list)
    }
}

# excute
run(config)
