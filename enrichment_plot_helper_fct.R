#library(gridExtra)
library(parallel)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(colorspace)
library(cowplot)

DCOLORPAL <- "RdBu"
#DCOLORPAL <- "Tropic"
SCOLORPAL <- "Plasma"

#source("~/Scripts/work/tools/map2color.r")
#source("~/Scripts/work/tools/iwanthue.r")



############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

tableGrob2 <- function(d, p = NULL) {
  # has_package("gridExtra")
  #d <- d[order(rownames(d)),]
  tp <- gridExtra::tableGrob(d)
  if (is.null(p)) {
    return(tp)
  }
  
  # Fix bug: The 'group' order of lines and dots/path is different
  p_data <- ggplot_build(p)$data[[1]]
  # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  p_data <- p_data[order(p_data[["group"]]), ]
  pcol <- unique(p_data[["colour"]])
  ## This is fine too
  ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
  j <- which(tp$layout$name == "rowhead-fg")
  
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
  }
  return(tp)
}

mygseaplot_multi <- function(xList, geneSetID, title = "", color = "green", base_size = 11, 
                             rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                             ES_geom = "line") 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(xList) == 1) {
    gsdata <- mygsInfo(xList[[1]], geneSetID)
    gsdata$COMP = names(xList)[1]
  }
  else {
    gsdata <- lapply(xList, mygsInfo, geneSetID = geneSetID)
    gsdata <- lapply(1:length(gsdata), function(j){
      gs <- gsdata[[j]]
      gs$COMP <- names(gsdata)[j]
      return(gs)
    })
    gsdata <- do.call(rbind, gsdata)
  }
  gsdata$COMP <- factor(gsdata$COMP, levels = names(xList))
  
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~COMP), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~COMP), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + geom_hline(yintercept=0, linetype="dashed", color = "black")
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  
  p.res <- p.res + geom_hline(yintercept=0, linetype="dashed", color = "black")
  
  i <- 0
  for (term in unique(gsdata$COMP)) {
    idx <- which(gsdata$ymin != 0 & gsdata$COMP == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~COMP)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) +
    theme(legend.position = "none", plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Rank in Ordered Dataset")
  
  
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
  if (length(color) == length(xList)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")# change color here
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    
    pdList <- lapply(1:length(xList), function(j){
      x <- xList[[j]]
      pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
      rownames(pd) <- names(xList)[j]
      pd <- pd[, -1]
      pd <- signif(pd, 2)
      return(pd)
    })
    pd <- do.call(rbind, pdList)
    
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                                                                                                                             0.95), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                    0.75), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                           0.9))
  }
  plotlist <- list(p.res, p2)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  aplot::plot_list(gglist = plotlist, ncol = 1, heights = rel_heights)
}

mygseaplot2 <- function(x, geneSetID, title = "", color = "green", base_size = 11, 
                        rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                        ES_geom = "line") 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- mygsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, mygsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + geom_hline(yintercept=0, linetype="dashed", color = "black")
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
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
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
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "NES", "pvalue", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    pd <- signif(pd, 2)
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.15), xmax = quantile(p.res$data$x, 
                                                                                                                              0.85), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                     0.75), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                            0.9))
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  aplot::plot_list(gglist = plotlist, ncol = 1, heights = rel_heights)
}

calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht, ht_gap = unit(1, "cm"), heatmap_legend_side = "top", annotation_legend_side = "top")
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  w = w + 1
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  h = h + 1
  dev.off()
  
  c(w, h)
}

getHeatmapMatrix <- function(gsea_list, TermIdx = 1, pvIdx = 4, qvIdx = NULL){
  # TermIdx: column index of the gene-set terms
  # pvIdx: column index of the pvalue or adjusted pvalue
  
  trms <- sort(unique(unlist(lapply(gsea_list, function(i) i[,TermIdx]))))
  m.pv <- matrix(1, nrow = length(trms), ncol = length(gsea_list))
  m.qv <- matrix(1, nrow = length(trms), ncol = length(gsea_list))
  
  for(i in 1:length(gsea_list))
  {
    gs <- gsea_list[[i]]
    idx <- match(unlist(gs[,TermIdx]), trms)
    m.pv[idx, i] <- as.numeric(gs[, pvIdx])
    
    if(!is.null(qvIdx)){
      m.qv[idx, i] <- as.numeric(gs[, qvIdx])
    }
  }
  rownames(m.pv) <- trms
  colnames(m.pv) <- names(gsea_list)
  
  rownames(m.qv) <- trms
  colnames(m.qv) <- names(gsea_list)
  return(list(pv=m.pv, qv=m.qv))
}

getScoreTwoWays <- function(inputList){
  # inputList is a list of pv matrix and qv matrix,
  # both should have the same size, with same row and column order.
  
  pvM <- inputList[["pv"]]
  qvM <- inputList[["qv"]]
  
  name.df <- data.frame(CNAME = colnames(pvM),
                        SIGN = ifelse(grepl(".UP", colnames(pvM)), "UP", "DOWN"),
                        ONEWAY = gsub(".UP|.DOWN", "", colnames(pvM)))
  oneway.name <- unique(name.df$ONEWAY)
  
  pvM.oneway <- matrix(1, nrow = nrow(pvM), ncol = length(oneway.name))
  qvM.oneway <- matrix(1, nrow = nrow(pvM), ncol = length(oneway.name))
  
  for(ri in 1:nrow(pvM)){
    for(ci in 1:length(oneway.name)){
      idx.up <- which(name.df$ONEWAY == oneway.name[ci] & name.df$SIGN == "UP") 
      idx.down <- which(name.df$ONEWAY == oneway.name[ci] & name.df$SIGN == "DOWN") 
      
      if(length(idx.up) != 0 & length(idx.down) != 0){
        if(pvM[ri, idx.up] <= pvM[ri, idx.down]){
          pvM.oneway[ri, ci] <- -log10(pvM[ri, idx.up])
          qvM.oneway[ri, ci] <- -log10(qvM[ri, idx.up])
        }else{
          pvM.oneway[ri, ci] <- +log10(pvM[ri, idx.down])
          qvM.oneway[ri, ci] <- +log10(qvM[ri, idx.down])
        }
      }else if(length(idx.up) != 0){
        pvM.oneway[ri, ci] <- -log10(pvM[ri, idx.up])
        qvM.oneway[ri, ci] <- -log10(qvM[ri, idx.up])
      }else{
        pvM.oneway[ri, ci] <- +log10(pvM[ri, idx.down])
        qvM.oneway[ri, ci] <- +log10(qvM[ri, idx.down])
      }
      
    }
  }
  rownames(pvM.oneway) <- rownames(pvM)
  colnames(pvM.oneway) <- oneway.name
  
  rownames(qvM.oneway) <- rownames(qvM)
  colnames(qvM.oneway) <- oneway.name
  
  return(list(pv=pvM.oneway, qv=qvM.oneway))
}

getScoreOneWay <- function(inputList){
  # inputList is a list of pv matrix and qv matrix,
  # both should have the same size, with same row and column order.
  
  return(list(pv=-log10(inputList[["pv"]]), qv=-log10(inputList[["qv"]])))
}

filterGsMatSingle <- function(gsMatList, nb, kw = NULL, gs = NULL, limit = 5, pvCutoff = 1, useQV = FALSE){
  
  gsMat <- gsMatList$pv
  
  if(!is.null(kw)) gsMat <- gsMat[grepl(kw, rownames(gsMat), ignore.case = TRUE), , drop = FALSE]
  if(!is.null(gs)) gsMat <- gsMat[match(intersect(rownames(gsMat), gs), rownames(gsMat)), , drop = FALSE]
  
  # select top X gene-sets per column
  idxList <- lapply(1:ncol(gsMat), function(x){
    values <- abs(as.numeric(gsMat[, x]))
    idx <- order(-values)
    idx <- idx[values[idx]>=-log10(pvCutoff)]
    if(length(idx) > nb) idx <- idx[1:nb]
    return(idx)
  })
  
  gsMat <- gsMat[unique(unlist(idxList)), , drop = FALSE]
  
  if(useQV) gsMat <- gsMatList$qv[match(rownames(gsMat), rownames(gsMatList$qv)), , drop = FALSE]
  
  gsMat[gsMat > limit] <- limit # set the limit to 5
  gsMat[gsMat < -limit] <- -limit
  return(gsMat)
}

geneSet_heatmap <- function(gsMat, outFile, colClust = FALSE){
  
  rowClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  
  myLimit <- max(abs(gsMat), na.rm = TRUE)
  if(myLimit == 0) myLimit <- 0.1
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = DCOLORPAL, reverse = TRUE)
  
  mywidth <- 0.75*ncol(gsMat)
  
  ifelse(nrow(gsMat) > 20, myheight <- 0.8*nrow(gsMat), myheight <- 0.8*nrow(gsMat))
  myheight <- max(c(myheight, 15))
  if(nrow(gsMat)<10) myheight <- 15
  
  gsMat.signif <- gsMat
  gsMat.signif[abs(gsMat) > -log10(0.05)] <- "*"
  gsMat.signif[gsMat.signif != "*"] <- ""
  
  h1 <- ComplexHeatmap::Heatmap(gsMat, name = "enrichment", col = col_fun, na_col = "black",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%s", gsMat.signif[i, j]), x, y, gp = gpar(fontsize = 14))},
                cluster_columns = colClust, cluster_rows = rowClust, 
                show_row_names = TRUE,
                column_title_rot = 90,
                row_names_max_width = unit(15, "cm"),
                column_names_max_height =  unit(15, "cm"),
                #row_names_gp = gpar(fontsize = 4),
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  pdfSize <- calc_ht_size(h1)
  pdf(outFile, width = pdfSize[1], height = pdfSize[2])
  draw(h1, ht_gap = unit(1, "cm"), heatmap_legend_side = "top", annotation_legend_side = "top")
  dev.off()
}

geneSet_heatmap_colAnn <- function(gsMat, outFile, colClust = FALSE, colAnn){
  
  colAnn <- colAnn[match(colnames(gsMat), rownames(colAnn)), , drop = FALSE]
  
  rowClust <- ifelse(nrow(gsMat) > 1, TRUE, FALSE)
  
  myLimit <- max(abs(gsMat), na.rm = TRUE)
  if(myLimit == 0) myLimit <- 0.1
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = DCOLORPAL, reverse = TRUE)
  
  mywidth <- 0.75*ncol(gsMat)
  
  ifelse(nrow(gsMat) > 20, myheight <- 0.8*nrow(gsMat), myheight <- 0.8*nrow(gsMat))
  myheight <- max(c(myheight, 15))
  if(nrow(gsMat)<10) myheight <- 15
  
  gsMat.signif <- gsMat
  gsMat.signif[abs(gsMat) > -log10(0.05)] <- "*"
  gsMat.signif[gsMat.signif != "*"] <- ""
  
  colCol <- lapply(colAnn, getColor)
  
  h1 <- Heatmap(gsMat, name = "enrichment", col = col_fun, na_col = "black",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%s", gsMat.signif[i, j]), x, y, gp = gpar(fontsize = 14))},
                cluster_columns = colClust, cluster_rows = rowClust, 
                show_row_names = TRUE,
                column_title_rot = 90,
                row_names_max_width = unit(15, "cm"),
                column_names_max_height =  unit(15, "cm"),
                
                column_split = colAnn[,1], cluster_column_slices = FALSE, column_gap = unit(5, "mm"),
                #row_split = myrowsplit, cluster_row_slices = FALSE, row_gap = unit(2, "mm"), row_title_rot = 0,
                top_annotation = HeatmapAnnotation(df = colAnn, col = colCol,
                                                   annotation_legend_param = list(direction = "horizontal")),
                
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  pdfSize <- calc_ht_size(h1)
  pdf(outFile, width = pdfSize[1], height = pdfSize[2])
  draw(h1, ht_gap = unit(1, "cm"), heatmap_legend_side = "top", annotation_legend_side = "top")
  dev.off()
}

heatmapPipe <- function(fh_list, outFile, nb = 5, TermIdx = 1, pvIdx = 4, qvIdx = NULL, kw = NULL, gs = NULL, limit = 5, pvCutoff =1, colClust = FALSE){
  mList <- getHeatmapMatrix(fh_list, TermIdx = TermIdx, pvIdx = pvIdx, qvIdx = qvIdx)
  mList.score <- getScoreTwoWays(mList)
  
  useQV = FALSE
  if(!is.null(qvIdx)) useQV = TRUE
  gsMat <- filterGsMatSingle(mList.score, nb = nb, kw = kw, gs = gs, limit = limit, pvCutoff = pvCutoff, useQV = useQV)
  
  geneSet_heatmap(gsMat, outFile, colClust = colClust)
}

heatmapPipeOneWay <- function(fh_list, outFile, nb = 5, TermIdx = 1, pvIdx = 4, qvIdx = NULL, kw = NULL, gs = NULL, limit = 5, pvCutoff =1, colClust = FALSE){
  mList <- getHeatmapMatrix(fh_list, TermIdx = TermIdx, pvIdx = pvIdx, qvIdx = qvIdx)
  mList.score <- getScoreOneWay(mList)
  
  useQV = FALSE
  if(!is.null(qvIdx)) useQV = TRUE
  gsMat <- filterGsMatSingle(mList.score, nb = nb, kw = kw, gs = gs, limit = limit, pvCutoff = pvCutoff, useQV = useQV)
  
  geneSet_heatmap(gsMat, outFile, colClust = colClust)
}

heatmapPipe_colAnn <- function(fh_list, outFile, nb = 5, TermIdx = 1, pvIdx = 4, qvIdx = NULL, kw = NULL, gs = NULL, limit = 5, pvCutoff =1, colClust = FALSE, colAnn){
  mList <- getHeatmapMatrix(fh_list, TermIdx = TermIdx, pvIdx = pvIdx, qvIdx = qvIdx)
  mList.score <- getScoreTwoWays(mList)
  
  useQV = FALSE
  if(!is.null(qvIdx)) useQV = TRUE
  gsMat <- filterGsMatSingle(mList.score, nb = nb, kw = kw, gs = gs, limit = limit, pvCutoff = pvCutoff, useQV = useQV)
  
  geneSet_heatmap_colAnn(gsMat, outFile, colClust = colClust, colAnn = colAnn)
}

overallBarplot <- function(fhFiles, nb, outFile, termIdx = 1, nbIdx = 11, pvIdx = 4, qvIdx = NULL, kw = NULL, gs = NULL, scoreMax = NULL, qvreorder = FALSE){
  fhList.up <- lapply(fhFiles, read.xlsx, sheet = "UP")
  fhList.down <- lapply(fhFiles, read.xlsx, sheet = "DOWN")
  
  # UP
  fh.up <- do.call(rbind, fhList.up)
  fh.up$"pvalue" <- as.numeric(fh.up[, pvIdx])
  fh.up$qvalue <- numeric(nrow(fh.up))
  if(nrow(fh.up)!=0) fh.up$"qvalue" <- 1
  if(!is.null(qvIdx)) fh.up$"qvalue" <-as.numeric(fh.up[, qvIdx])
  fh.up$"NB" <- as.numeric(fh.up[, nbIdx])
  fh.up <- fh.up[order(fh.up$"pvalue"), ]
  colnames(fh.up)[termIdx] <- "Term"
  
  if(!is.null(kw)) fh.up <- fh.up[grepl(kw, fh.up$Term, ignore.case = TRUE), ]
  if(!is.null(gs)) fh.up <- fh.up[match(intersect(fh.up$Term, gs), fh.up$Term), ]
  
  if(nrow(fh.up) > nb){
    ggmat.up <- fh.up[1:nb, c("Term", "pvalue", "qvalue", "NB")]
  }else{
    ggmat.up <- fh.up[, c("Term", "pvalue", "qvalue", "NB")]
  }
  
  fh.down <- do.call(rbind, fhList.down)
  fh.down$"pvalue" <- as.numeric(fh.down[, pvIdx])
  fh.down$qvalue <- numeric(nrow(fh.down))
  if(nrow(fh.down)!=0) fh.down$"qvalue" <- 1
  if(!is.null(qvIdx)) fh.down$"qvalue" <-as.numeric(fh.down[, qvIdx])
  fh.down$"NB" <- as.numeric(fh.down[, nbIdx])
  fh.down <- fh.down[order(fh.down$"pvalue"), ]
  colnames(fh.down)[termIdx] <- "Term"
  
  if(!is.null(kw)) fh.down <- fh.down[grepl(kw, fh.down$Term, ignore.case = TRUE), ]
  if(!is.null(gs)) fh.down <- fh.down[match(intersect(fh.down$Term, gs), fh.down$Term), ]

  if(nrow(fh.down) > nb){
    ggmat.down <- fh.down[1:nb, c("Term", "pvalue", "qvalue", "NB")]
  }else{
    ggmat.down <- fh.down[, c("Term", "pvalue", "qvalue", "NB")]
  }
  
  # ADJUST nb
  nb <- min(c(nb, max(c(nrow(ggmat.up), nrow(ggmat.down)))))
  
  # SET Y-AXIS LIMIT (score = -log10 pvalue)
  if(is.null(scoreMax)){
    scoreMax <- max(ceiling(c(-log10(ggmat.up$"pvalue"), -log10(ggmat.down$"pvalue"))))
    if(!is.null(qvIdx)) scoreMax <- max(ceiling(c(-log10(ggmat.up$"qvalue"), -log10(ggmat.down$"qvalue"))))
    if(scoreMax < 2) scoreMax <- 2
  }
  
  # Add fake data to fill empty space
  if(nrow(ggmat.up) != 0 & nrow(ggmat.up) < nb){
    nb.missing <- nb - nrow(ggmat.up)
    ggmat.fake <- data.frame("Term" = paste0("NA", 1:nb.missing),
                             "pvalue" = 1,
                             "qvalue" = 1,
                             "NB" = min(ggmat.up$"NB"))
    ggmat.up <- rbind(ggmat.up, ggmat.fake)
  }
  
  ggmat.up$Score <- -log10(ggmat.up$"pvalue")
  if(!is.null(qvIdx)) ggmat.up$Score <- -log10(ggmat.up$"qvalue")
  ggmat.up$Score[ggmat.up$Score > scoreMax] <- scoreMax
  ggmat.up$Term <- gsub('(.{1,70})', '\\1\n', ggmat.up$Term)
  ggmat.up$Term <- factor(ggmat.up$Term, levels = rev(ggmat.up$Term[order(ggmat.up$"pvalue")]))
  if(qvreorder) ggmat.up$Term <- factor(ggmat.up$Term, levels = rev(ggmat.up$Term[order(ggmat.up$"qvalue")]))
  
  if(scoreMax < 5){
    mybreaks <- seq(0, scoreMax, by = 1)
  }else if(scoreMax < 10){
    mybreaks <- seq(0, scoreMax, by = 2)
  }else if(scoreMax < 30){
    mybreaks <- seq(0, scoreMax, by = 5)
  }else{
    mybreaks <- seq(0, scoreMax, by = 10)
  }
  
  p <- ggplot(data=ggmat.up, aes(x=Term, y=Score, fill = NB))
  p <- p +  geom_bar(stat="identity", color = "black")
  p <- p + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey20", linewidth = 1.25)
  p <- p + coord_flip()
  p <- p + scale_y_continuous(limits = c(0, scoreMax), breaks = mybreaks)
  p <- p + scale_fill_continuous_sequential(palette = "Reds")
  p <- p + theme_bw(base_size = 10)
  p <- p + theme(text=element_text(color="black"),axis.text=element_text(color="black"))
  p <- p + xlab("")
  p <- p + ggtitle("UP")
  if(!is.null(qvIdx)){
    p <- p + ylab("-log10 qvalue")
  }else{
    p <- p + ylab("-log10 pvalue")
  }
  
  p.up <- p
  #ggsave(plot = p, filename = paste0(outName, "_top", nb, "UP_overall_barplot.pdf"), width = 8, height = 10)
  
  # down
  
  # Add fake data to fill empty space
  if(nrow(ggmat.down) != 0 & nrow(ggmat.down) < nb){
    nb.missing <- nb - nrow(ggmat.down)
    ggmat.fake <- data.frame("Term" = paste0("NA", 1:nb.missing),
                             "pvalue" = 1,
                             "qvalue" = 1,
                             "NB" = min(ggmat.down$"NB"))
    ggmat.down <- rbind(ggmat.down, ggmat.fake)
  }
  
  ggmat.down$Score <- -log10(ggmat.down$"pvalue")
  if(!is.null(qvIdx)) ggmat.down$Score <- -log10(ggmat.down$"qvalue")
  ggmat.down$Score[ggmat.down$Score > scoreMax] <- scoreMax
  ggmat.down$Term <- gsub('(.{1,70})', '\\1\n', ggmat.down$Term)
  ggmat.down$Term <- factor(ggmat.down$Term, levels = rev(ggmat.down$Term[order(ggmat.down$"pvalue")]))
  if(qvreorder) ggmat.down$Term <- factor(ggmat.down$Term, levels = rev(ggmat.down$Term[order(ggmat.down$"qvalue")]))
  
  
  p <- ggplot(data=ggmat.down, aes(x=Term, y=Score, fill = NB))
  p <- p +  geom_bar(stat="identity", color = "black")
  p <- p + geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey20", linewidth = 1.25)
  p <- p + coord_flip()
  p <- p + scale_y_continuous(limits = c(0, scoreMax), breaks = mybreaks)
  p <- p + scale_fill_continuous_sequential(palette = "Blues")
  p <- p + theme_bw(base_size = 10)
  p <- p + theme(text=element_text(color="black"),axis.text=element_text(color="black"))
  p <- p + xlab("")
  p <- p + ggtitle("DOWN")
  if(!is.null(qvIdx)){
    p <- p + ylab("-log10 qvalue")
  }else{
    p <- p + ylab("-log10 pvalue")
  }
  
  p.down <- p
  #ggsave(plot = p, filename = paste0(outName, "_top", nb, "DOWN_overall_barplot.pdf"), width = 8, height = 10)
  
  #p <- arrangeGrob(p.up, p.down, ncol = 2) #generates g
  p <- plot_grid(p.up, p.down, labels = NULL)
  myheight <- 3
  if(nb > 10) myheight <- 2 + (nb*0.15)
  #ggsave(plot = p, filename = paste0(outName, "_top", nb, "_overall_barplot.pdf"), width = 8*2, height = myheight)
  ggsave(plot = p, filename = outFile, width = 8*2, height = myheight)
  
}

getNESmat <- function(fhList, termIdx = 1, nesIdx = 4){
  
  geneSets <- sort(unique(unlist(lapply(fhList, function(x) return(as.character(x[, termIdx]))))))
  
  nesMat <- matrix(0, nrow = length(geneSets), ncol = length(fhList))
  
  for(yi in 1:ncol(nesMat)){
    fh.current <- fhList[[yi]]
    nesMat[, yi] <- fh.current[match(geneSets, as.character(fh.current[, termIdx])), nesIdx]
  }
  nesMat[is.na(nesMat)] <- 0
  rownames(nesMat) <- geneSets
  colnames(nesMat) <- names(fhList)
  
  return(data.frame(nesMat))
}


