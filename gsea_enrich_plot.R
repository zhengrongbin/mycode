
library(fgsea)
library(ggplot2)

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
    for (geneSetId in geneSetID){
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

enrich_plot = function(geneList, geneSet, geneSetID, pre_ranked = FALSE,
    title = "", color = "green", base_size = 11 ,
    rel_heights = c(1.5, 0.5, 1), subplots = 1:3,
    pvalue_table = FALSE, ES_geom = "line")
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
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
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
    col <- c(rev(RColorBrewer::brewer.pal(5, "Blues")), RColorBrewer::brewer.pal(5, "Reds"))
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
    
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
    axis.ticks.x = element_line(), axis.text.x = element_text())
    if (length(subplots) == 1) {
        return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
            r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
    }
    if (length(rel_heights) > length(subplots)) {
        rel_heights <- rel_heights[subplots]
    }
    cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}
