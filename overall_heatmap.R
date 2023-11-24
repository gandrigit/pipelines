
library(org.Hs.eg.db)
library(org.Hs.eg.db)
library(openxlsx)
library(ComplexHeatmap)
library(viridis)
library(GeneAnswers)
library(circlize)
library(measurements)
library(limma)
library(msigdbr)
library(circlize)
library(parallel)

# EXAMPLE OF COLOR PALETTES
if(FALSE){
  # COLOR PALETTES
  hcl.pals("diverging")# "Blue-Red 2"
  hcl.pals("divergingx")# "PuOr" "RdBu" "Tropic"
  hcl.pals("sequential")# "Plasma"
}

# some other sequential ones from archR
horizon = c("1"='#000075',"4"='#2E00FF', "6"='#9408F7', "10"='#C729D6', "8"='#FA4AB5', "3"='#FF6A95', "7"='#FF8B74', "5"='#FFAC53', "9"='#FFCD32', "2"='#FFFF60')
horizonExtra =c("1"="#000436","4"="#021EA9","6"="#1632FB","8"="#6E34FC","3"="#C732D5","9"="#FD619D","7"="#FF9965","5"="#FFD32B","2"="#FFFC5A")
solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')
comet = c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black")
fireworks2 = c("5"="black", "2"="#2488F0","4"="#7F3F98","3"="#E22929","1"="#FCB31A")
greyMagma = c("2"="#E6E7E8", "4"="#FB8861FF", "5"="#B63679FF", "3"="#51127CFF", "1"="#000004FF")


# Consensus Path DB
#load("~/Research/database/consensus/Consensus_unique08.mouse.RData") # load cons2 annotation list
# Misgdb 
#load("~/Research/database/msigdb_mouse/msigdb_v7.0.mouse.RData")
#dbList[["Consensus"]] <- consdb

# BUILD DB
m_df <- msigdbr(species = "Homo sapiens")
m_df$gs_subcat <- gsub("\\:", "_", m_df$gs_subcat)
m_df$gs_subcat <- paste0(m_df$gs_cat, "_", m_df$gs_subcat)
m_df$gs_subcat <- gsub("_$", "", m_df$gs_subcat)

subcat.unique <- sort(unique(m_df$gs_subcat))

dbList <- lapply(subcat.unique, function(i){
  m_df.sub <- m_df[m_df$gs_subcat == i, ]
  terms <- unique(m_df.sub$gs_name)
  db <- mclapply(terms, function(j){
    return(m_df.sub$entrez_gene[m_df.sub$gs_name == j])
  }, mc.cores = 8)
  names(db) <- terms
  return(db)
})
names(dbList) <- subcat.unique

#dir.create("~/Research/database/misgdbr/")
#save(dbList, file = "~/Research/database/misgdbr/dbList_human_7.2.1.RData")

# ADD CONSENSUS
#load(file.path("~/Research/database/consensus_path_db/release34/consensus34.RData"))
#dbList[["Consensus"]] <- cons34

source("~/Scripts/work/tools/iwanthue.r")

############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################



my.read_xlsx <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile)
  names(mList) <- mysheets
  return(mList)
}

my.read_xlsxNoHeader <- function(inFile)
{
  mysheets <- getSheetNames(inFile)
  mList <- lapply(mysheets, read.xlsx, xlsxFile = inFile, colNames = FALSE)
  names(mList) <- mysheets
  return(mList)
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unlist(lapply(symbol, function(i) return(i[1])))
  return(symbol)
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

getUPandDOWN <- function(degMat, nb){
  degMat.up <- degMat[degMat$fc > 0, ]
  degMat.up <- degMat.up[order(degMat.up$pv), ]
  
  degMat.down <- degMat[degMat$fc < 0, ]
  degMat.down <- degMat.down[order(degMat.down$pv), ]
  
  up <- degMat.up$entrez
  if(length(up) > nb) up <- up[1:nb]
  
  down <- degMat.down$entrez
  if(length(down) > nb) down <- down[1:nb]
  
  return(c(up, down))
}

robust_dist = function(x, y) {
  qx = quantile(x, c(0.1, 0.9))
  qy = quantile(y, c(0.1, 0.9))
  l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
  x = x[l]
  y = y[l]
  sqrt(sum((x - y)^2))
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


# INTENSITY AND FOLD CHANGE
overall_heatmap <- function(expr, annot, degList, mygenes, pv = 0.05, fc = 0, nbMAX = NULL, sortIdx = NULL, onlyDEG = FALSE, legend.pos = "bottom", outFile){
  
  if(onlyDEG){
    degList.deg <- lapply(degList, function(k) k$entrez[k$qv < pv & abs(k$fc) > fc])
    mygenes <- intersect(mygenes, unlist(degList.deg))
  }
  
  if(length(mygenes) == 0){
    print(paste0("no genes in ", outFile))
    return(NA)
  }
  
  # H2: fc
  degList.sort <- lapply(degList, function(k) k[order(k$entrez), ])
  fcList <- lapply(degList.sort, function(k) k$fc)
  fcMat <- do.call(cbind, fcList)
  rownames(fcMat) <- degList.sort[[1]]$entrez
  colnames(fcMat) <- names(degList)
  
  fcMat <- fcMat[match(mygenes, rownames(fcMat)), , drop = FALSE]
  rownames(fcMat) <- degList[[1]]$symbol[match(rownames(fcMat), degList[[1]]$entrez)]
  
  qvList <- lapply(degList.sort, function(k) k$qv)
  qvMat <- do.call(cbind, qvList)
  rownames(qvMat) <- degList.sort[[1]]$entrez
  colnames(qvMat) <- names(degList)
  
  qvMat <- qvMat[match(mygenes, rownames(qvMat)), , drop = FALSE]
  qvMat.signif <- qvMat
  qvMat.signif[qvMat < pv & abs(fcMat) > fc] <- "*"
  qvMat.signif[qvMat.signif != "*"] <- ""
  
  # SELECT THE TOP XXX GENES BASED ON MINIMUM PVALUE AND THEN MAX ABSOLUTE LOGFC
  if(!is.null(nbMAX) & nrow(fcMat) > nbMAX){
    qvMat.min <- apply(qvMat, 1, min)
    fcMat.abs.max <- apply(abs(fcMat), 1, max)
    
    topOrder <- order(qvMat.min, -fcMat.abs.max)[1:nbMAX]
    fcMat <- fcMat[topOrder, , drop = FALSE]
    qvMat <- qvMat[topOrder, , drop = FALSE]
    qvMat.signif <- qvMat.signif[topOrder, , drop = FALSE]
    mygenes <- mygenes[topOrder]
  }
  
  # RE-ORDER
  if(!is.null(sortIdx)){
    doClust <- FALSE
    newOrder <- order(-as.numeric(fcMat[, sortIdx]))
    fcMat <- fcMat[newOrder, , drop = FALSE]
    qvMat <- qvMat[newOrder, , drop = FALSE]
    qvMat.signif <- qvMat.signif[newOrder, , drop = FALSE]
  }else{
    doClust <- TRUE
  }
  
  
  myLimit <- max(abs(fcMat), na.rm = TRUE)
  myLimit <- max(c(1, myLimit))
  #col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("skyblue", "white", "tomato"))
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = "PuOr", reverse = TRUE)
  
  mywidth <- 0.75*ncol(fcMat)
  
  ifelse(nrow(fcMat) > 20, myheight <- 0.8*nrow(fcMat), myheight <- 1.5*nrow(fcMat))
  myheight <- max(c(myheight, 15))
  
  h2 <- Heatmap(fcMat, name = "FC", col = col_fun, na_col = "black",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%s", qvMat.signif[i, j]), x, y, gp = gpar(fontsize = 14))},
                cluster_columns = FALSE, cluster_rows = doClust, 
                column_names_side = "top",
                show_row_names = TRUE,
                row_names_max_width = unit(6, "cm"),
                column_names_max_height =  unit(6, "cm"),
                column_title_rot = 90,
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  # H1: zscore expression
  intMat <- expr[mygenes, , drop = FALSE]
  rownames(intMat) <- degList[[1]]$symbol[match(rownames(intMat), degList[[1]]$entrez)]
  
  # RE-ORDER
  if(!is.null(sortIdx)){
    doClust <- FALSE
    intMat <- intMat[newOrder, , drop = FALSE]
  }else{
    doClust <- TRUE
  }
  
  # Column annotation
  phAnn.col <- annot[match(colnames(intMat), annot$SAMPLE), c("GROUP"), drop = FALSE] # change this!
  rownames(phAnn.col) <- colnames(intMat)
  
  colCol <- lapply(phAnn.col, getColor)
  
  myLimit <- max(abs(intMat), na.rm = TRUE)
  #col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("royalblue", "snow", "firebrick1"))
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = "RdBu", reverse = TRUE)
  
  mywidth <- 0.75*ncol(intMat)
  
  h1 <- Heatmap(intMat, name = "Zscore", col = col_fun, na_col = "black",
                cluster_columns = TRUE, cluster_rows = doClust, 
                show_row_names = TRUE,
                row_names_max_width = unit(6, "cm"),
                column_names_max_height =  unit(6, "cm"),
                column_title_rot = 90,
                column_split = phAnn.col[,1], cluster_column_slices = FALSE, column_gap = unit(5, "mm"),
                #row_split = myrowsplit, cluster_row_slices = FALSE, row_gap = unit(2, "mm"), row_title_rot = 0,
                top_annotation = HeatmapAnnotation(df = phAnn.col, col = colCol,
                                                   annotation_legend_param = list(direction = "horizontal")),
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  
  ht_list <- h1 + h2
  #draw(ht_list, ht_gap = unit(1, "cm"))
  
  # save pdf
  #pdfHeight <- conv_unit(as.numeric(ComplexHeatmap:::height(draw(ht_list, ht_gap = unit(1, "cm")))), "mm", "inch")+5
  #pdfHeight <- max(c(pdfHeight, 8))
  #pdfWidth <- conv_unit(as.numeric(ComplexHeatmap:::width(draw(ht_list, ht_gap = unit(1, "cm")))), "mm", "inch")+10
  
  pdfSize <- calc_ht_size(ht_list)
  
  pdf(outFile, width = pdfSize[1], height = pdfSize[2])
  draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = legend.pos, annotation_legend_side = legend.pos)
  dev.off()
}

overall_heatmap_NOLABEL <- function(expr, annot, degList, mygenes, nbMAX = NULL, sortIdx = NULL, onlyDEG = FALSE, legend.pos = "bottom", outFile){
  
  if(onlyDEG){
    degList.deg <- lapply(degList, function(k) k$entrez[k$qv < pv & abs(k$fc) > fc])
    mygenes <- intersect(mygenes, unlist(degList.deg))
  }
  
  if(length(mygenes) == 0){
    print(paste0("no genes in ", outFile))
    return(NA)
  }
  
  # H2: fc
  degList.sort <- lapply(degList, function(k) k[order(k$entrez), ])
  fcList <- lapply(degList.sort, function(k) k$fc)
  fcMat <- do.call(cbind, fcList)
  rownames(fcMat) <- degList.sort[[1]]$entrez
  colnames(fcMat) <- names(degList)
  
  fcMat <- fcMat[match(mygenes, rownames(fcMat)), , drop = FALSE]
  rownames(fcMat) <- degList[[1]]$symbol[match(rownames(fcMat), degList[[1]]$entrez)]
  
  qvList <- lapply(degList.sort, function(k) k$qv)
  qvMat <- do.call(cbind, qvList)
  rownames(qvMat) <- degList.sort[[1]]$entrez
  colnames(qvMat) <- names(degList)
  
  qvMat <- qvMat[match(mygenes, rownames(qvMat)), , drop = FALSE]
  
  # SELECT THE TOP XXX GENES BASED ON MINIMUM PVALUE AND THEN MAX ABSOLUTE LOGFC
  if(!is.null(nbMAX) & nrow(fcMat) > nbMAX){
    qvMat.min <- apply(qvMat, 1, min)
    fcMat.abs.max <- apply(abs(fcMat), 1, max)
    
    topOrder <- order(qvMat.min, -fcMat.abs.max)[1:nbMAX]
    fcMat <- fcMat[topOrder, , drop = FALSE]
    qvMat <- qvMat[topOrder, , drop = FALSE]
    mygenes <- mygenes[topOrder]
  }
  
  # RE-ORDER
  if(!is.null(sortIdx)){
    doClust <- FALSE
    newOrder <- order(-as.numeric(fcMat[, sortIdx]))
    fcMat <- fcMat[newOrder, , drop = FALSE]
    qvMat <- qvMat[newOrder, , drop = FALSE]
  }else{
    doClust <- TRUE
  }
  
  
  myLimit <- max(abs(fcMat), na.rm = TRUE)
  myLimit <- max(c(1, myLimit))
  #col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("skyblue", "white", "tomato"))
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = "PuOr", reverse = TRUE)
  
  mywidth <- 0.75*ncol(fcMat)
  
  ifelse(nrow(fcMat) > 20, myheight <- 0.8*nrow(fcMat), myheight <- 1.5*nrow(fcMat))
  myheight <- min(c(myheight, 10))
  
  h2 <- Heatmap(fcMat, name = "FC", col = col_fun, na_col = "black",
                cluster_columns = FALSE, cluster_rows = doClust, 
                column_names_side = "top",
                show_row_names = FALSE,
                column_title_rot = 90,
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  # H1: zscore expression
  intMat <- expr[mygenes, , drop = FALSE]
  rownames(intMat) <- degList[[1]]$symbol[match(rownames(intMat), degList[[1]]$entrez)]
  
  # RE-ORDER
  if(!is.null(sortIdx)){
    doClust <- FALSE
    intMat <- intMat[newOrder, , drop = FALSE]
  }else{
    doClust <- TRUE
  }
  
  # Column annotation
  phAnn.col <- annot[match(colnames(intMat), annot$SAMPLE), c("GROUP"), drop = FALSE] # change this!
  rownames(phAnn.col) <- colnames(intMat)
  
  colCol <- lapply(phAnn.col, getColor)
  
  myLimit <- max(abs(intMat), na.rm = TRUE)
  #col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("royalblue", "snow", "firebrick1"))
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = "RdBu", reverse = TRUE)
  
  mywidth <- 0.75*ncol(intMat)
  
  h1 <- Heatmap(intMat, name = "Zscore", col = col_fun, na_col = "black",
                cluster_columns = TRUE, cluster_rows = doClust, 
                show_row_names = FALSE,
                column_title_rot = 90,
                column_split = phAnn.col[,1], cluster_column_slices = FALSE, column_gap = unit(5, "mm"),
                #row_split = myrowsplit, cluster_row_slices = FALSE, row_gap = unit(2, "mm"), row_title_rot = 0,
                top_annotation = HeatmapAnnotation(df = phAnn.col, col = colCol,
                                                   annotation_legend_param = list(direction = "horizontal")),
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  
  
  ht_list <- h1 + h2
  
  # save pdf
  pdfSize <- calc_ht_size(ht_list)
  
  pdf(outFile, width = pdfSize[1], height = pdfSize[2])
  draw(ht_list, ht_gap = unit(1, "cm"), heatmap_legend_side = legend.pos, annotation_legend_side = legend.pos)
  dev.off()
  
  
}

# FOLD CHANGE ONLY
fc_heatmap <- function(annot, degList, mygenes, pv = 0.05, fc = 0, nbMAX = NULL, sortIdx = NULL, onlyDEG = FALSE, legend.pos = "bottom", outFile){
  
  if(onlyDEG){
    degList.deg <- lapply(degList, function(k) k$entrez[k$qv < pv & abs(k$fc) > fc])
    mygenes <- intersect(mygenes, unlist(degList.deg))
  }
  
  if(length(mygenes) == 0){
    print(paste0("no genes in ", outFile))
    return(NA)
  }
  
  # H2: fc
  degList.sort <- lapply(degList, function(k) k[order(k$entrez), ])
  fcList <- lapply(degList.sort, function(k) k$fc)
  fcMat <- do.call(cbind, fcList)
  rownames(fcMat) <- degList.sort[[1]]$entrez
  colnames(fcMat) <- names(degList)
  
  fcMat <- fcMat[match(mygenes, rownames(fcMat)), , drop = FALSE]
  rownames(fcMat) <- degList[[1]]$symbol[match(rownames(fcMat), degList[[1]]$entrez)]
  
  qvList <- lapply(degList.sort, function(k) k$qv)
  qvMat <- do.call(cbind, qvList)
  rownames(qvMat) <- degList.sort[[1]]$entrez
  colnames(qvMat) <- names(degList)
  
  qvMat <- qvMat[match(mygenes, rownames(qvMat)), , drop = FALSE]
  qvMat.signif <- qvMat
  qvMat.signif[qvMat < pv & abs(fcMat) > fc] <- "*"
  qvMat.signif[qvMat.signif != "*"] <- ""
  
  # SELECT THE TOP XXX GENES BASED ON MINIMUM PVALUE AND THEN MAX ABSOLUTE LOGFC
  if(!is.null(nbMAX) & nrow(fcMat) > nbMAX){
    qvMat.min <- apply(qvMat, 1, min)
    fcMat.abs.max <- apply(abs(fcMat), 1, max)
    
    topOrder <- order(qvMat.min, -fcMat.abs.max)[1:nbMAX]
    fcMat <- fcMat[topOrder, , drop = FALSE]
    qvMat <- qvMat[topOrder, , drop = FALSE]
    qvMat.signif <- qvMat.signif[topOrder, , drop = FALSE]
    mygenes <- mygenes[topOrder]
  }
  
  # RE-ORDER
  if(!is.null(sortIdx)){
    doClust <- FALSE
    newOrder <- order(-as.numeric(fcMat[, sortIdx]))
    fcMat <- fcMat[newOrder, , drop = FALSE]
    qvMat <- qvMat[newOrder, , drop = FALSE]
    qvMat.signif <- qvMat.signif[newOrder, , drop = FALSE]
  }else{
    doClust <- TRUE
  }
  
  
  myLimit <- max(abs(fcMat), na.rm = TRUE)
  myLimit <- max(c(1, myLimit))
  #col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("skyblue", "white", "tomato"))
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = "PuOr", reverse = TRUE)
  
  mywidth <- 0.75*ncol(fcMat)
  
  ifelse(nrow(fcMat) > 20, myheight <- 0.8*nrow(fcMat), myheight <- 1.5*nrow(fcMat))
  myheight <- max(c(myheight, 15))
  
  h2 <- Heatmap(fcMat, name = "FC", col = col_fun, na_col = "black",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%s", qvMat.signif[i, j]), x, y, gp = gpar(fontsize = 14))},
                cluster_columns = FALSE, cluster_rows = doClust, 
                column_names_side = "top",
                show_row_names = TRUE,
                column_title_rot = 90,
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  # save pdf
  pdfSize <- calc_ht_size(h2)
  
  pdf(outFile, width = pdfSize[1], height = pdfSize[2])
  draw(h2, ht_gap = unit(1, "cm"), heatmap_legend_side = legend.pos, annotation_legend_side = legend.pos)
  dev.off()
  
  
}

fc_heatmap_NOLABEL <- function(annot, degList, mygenes, nbMAX = NULL, sortIdx = NULL, onlyDEG = FALSE, legend.pos = "bottom", outFile){
  
  if(onlyDEG){
    degList.deg <- lapply(degList, function(k) k$entrez[k$qv < pv & abs(k$fc) > fc])
    mygenes <- intersect(mygenes, unlist(degList.deg))
  }
  
  if(length(mygenes) == 0){
    print(paste0("no genes in ", outFile))
    return(NA)
  }
  
  # H2: fc
  degList.sort <- lapply(degList, function(k) k[order(k$entrez), ])
  fcList <- lapply(degList.sort, function(k) k$fc)
  fcMat <- do.call(cbind, fcList)
  rownames(fcMat) <- degList.sort[[1]]$entrez
  colnames(fcMat) <- names(degList)
  
  fcMat <- fcMat[match(mygenes, rownames(fcMat)), , drop = FALSE]
  rownames(fcMat) <- degList[[1]]$symbol[match(rownames(fcMat), degList[[1]]$entrez)]
  
  qvList <- lapply(degList.sort, function(k) k$qv)
  qvMat <- do.call(cbind, qvList)
  rownames(qvMat) <- degList.sort[[1]]$entrez
  colnames(qvMat) <- names(degList)
  
  qvMat <- qvMat[match(mygenes, rownames(qvMat)), , drop = FALSE]
  
  # SELECT THE TOP XXX GENES BASED ON MINIMUM PVALUE AND THEN MAX ABSOLUTE LOGFC
  if(!is.null(nbMAX) & nrow(fcMat) > nbMAX){
    qvMat.min <- apply(qvMat, 1, min)
    fcMat.abs.max <- apply(abs(fcMat), 1, max)
    
    topOrder <- order(qvMat.min, -fcMat.abs.max)[1:nbMAX]
    fcMat <- fcMat[topOrder, , drop = FALSE]
    qvMat <- qvMat[topOrder, , drop = FALSE]
    mygenes <- mygenes[topOrder]
  }
  
  # RE-ORDER
  if(!is.null(sortIdx)){
    doClust <- FALSE
    newOrder <- order(-as.numeric(fcMat[, sortIdx]))
    fcMat <- fcMat[newOrder, , drop = FALSE]
    qvMat <- qvMat[newOrder, , drop = FALSE]
  }else{
    doClust <- TRUE
  }
  
  
  myLimit <- max(abs(fcMat), na.rm = TRUE)
  myLimit <- max(c(1, myLimit))
  #col_fun <- colorRamp2(c(-myLimit, 0, myLimit), c("skyblue", "white", "tomato"))
  col_fun <- colorRamp2(c(-myLimit, 0, myLimit), hcl_palette = "PuOr", reverse = TRUE)
  
  mywidth <- 0.75*ncol(fcMat)
  
  ifelse(nrow(fcMat) > 20, myheight <- 0.8*nrow(fcMat), myheight <- 1.5*nrow(fcMat))
  myheight <- min(c(myheight, 10))
  
  h2 <- Heatmap(fcMat, name = "FC", col = col_fun, na_col = "black",
                cluster_columns = FALSE, cluster_rows = doClust, 
                column_names_side = "top",
                show_row_names = FALSE,
                column_title_rot = 90,
                width = unit(mywidth, "cm"), heatmap_height = unit(myheight , "cm"),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  
  
  # save pdf
  pdfSize <- calc_ht_size(h2)
  
  pdf(outFile, width = pdfSize[1], height = pdfSize[2])
  draw(h2, ht_gap = unit(1, "cm"), heatmap_legend_side = legend.pos, annotation_legend_side = legend.pos)
  dev.off()
  
  
}



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

###########################################################
# PARAMETERS
docDir <- file.path("~/Research/XXX/doc")

countDir <- file.path("~/Research/XXX/count/")

degDir <- file.path("~/Research/XXX/limma_voom/")

heatmapDir.main <- file.path("~/Research/XXX/overallHeatmaps/")
dir.create(heatmapDir.main, recursive = TRUE, showWarnings = FALSE)


###########################################################
# LOAD EXPRESSION
setwd(countDir)
#expr <- read.delim("rlogTransformation.txt")# DESeq2
expr <- read.delim("logCPM.txt")# DESeq2


# ANNOTATIONS
setwd(docDir)
ann.sample <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann.sample <- ann.sample[match(colnames(expr), ann.sample$SAMPLE), ]


# ADD GENE META DATA
ensID <- rownames(expr)
entrezID <- ensembl2entrez(ensID)

expr <- expr[!is.na(entrezID), ]
ensID <- ensID[!is.na(entrezID)]
entrezID <- entrezID[!is.na(entrezID)]

#  ENSEMBL TO ENTREZ (IQR)
expr.IQR <- apply(expr, 1, IQR)
expr <- lapply(unique(entrezID), function(i){
  ensID.current <- ensID[entrezID == i]
  if(length(ensID.current) == 1) expr[ensID.current,] else{
    expr[names(which.max(expr.IQR[ensID.current])), ]
  }
})
expr <- do.call(rbind, expr)
rownames(expr) <- unique(entrezID)


# SET GROUPS
ann.sample$GROUP <- factor(ann.sample$GROUP, levels = c("G1", "G2"
))

# RE-ORDER SAMPLES
new_order <- order(ann.sample$GROUP)
ann.sample <- ann.sample[new_order, ]
expr <- expr[, new_order]

# REMOVE GROUP
#toRemove <- ann.sample$GROUP == "Con_HLH"
#ann.sample <- ann.sample[!toRemove, ]
#expr <- expr[, !toRemove]


###########################################################
# LOAD LIMMA  
setwd(degDir)
degFiles <- list.files(pattern = "_limma.xlsx")
names(degFiles) <- gsub("_limma.xlsx", "", degFiles)
degList <- lapply(degFiles, read.xlsx, sheet = 1)

degList <- degList[c("X-Y", "Y-Z", "...")]

degList <- lapply(degList, function(i){
  i <- i[!is.na(i$entrez), ]
  i <- i[!duplicated(i$entrez), ]
  i <- i[!is.na(i$symbol), ]
  i <- i[!duplicated(i$symbol), ]
  
  df <- data.frame(entrez = i$entrez,
                   symbol = i$symbol,
                   fc = i$logFC,
                   pv = i$P.Value,
                   qv = i$adj.P.Val)
  return(df)
})


# INTERSECT expr and deseqList
common <- intersect(rownames(expr), degList[[1]]$entrez)
expr <- expr[common, ]
degList <- lapply(degList, function(i) i[i$entrez %in% common, ])


# SCALE DATA
expr.scaled <- t(scale(t(expr)))

if(FALSE){
  gene.meta <- data.frame(entrez = rownames(expr),
                          symbol = degList[[1]]$symbol[match(rownames(expr), degList[[1]]$entrez)])
  
  toxlsx <- list(log2CPM = data.frame(cbind(gene.meta, expr)),
                 zscore = data.frame(cbind(gene.meta, t(scale(t(expr))))),
                 annotation = ann.sample)
  write.xlsx(toxlsx, file.path(countDir, "normalized_annotated.xlsx"), rowNames = FALSE)
}

pvCutoff <- 0.05
fcCutoff <- 0

###########################################################
# LOOP OVER DBLIST

dbList.sub <- dbList[c("H")]

topList <- lapply(degList, getUPandDOWN, 20)
dbList.sub <- list("DEG" = topList)

if(FALSE){
  # OR BUILD IT FROM KW
  togrep <- c("GSE30962")
  dbList.sub <- lapply(togrep, function(j){
    sub <- lapply(dbList, function(i){
      i[base::grepl(j, names(i), ignore.case = TRUE)]
    })
    sub <- sub[sapply(sub, length) != 0]
    sub <- do.call(c, sub)
  })
  names(dbList.sub) <- togrep
}

#
lapply(1:length(dbList.sub), function(i){
  gsList <- dbList.sub[[i]]
  gsDirName <- names(dbList.sub)[i]
  heatmapDir <- file.path(heatmapDir.main, gsDirName)
  dir.create(heatmapDir, recursive = TRUE)
  
  # get genes-of-interest
  lapply(1:length(gsList), function(j){
    gsName <- gsub("\\/", "__", names(gsList)[j])
    gsName <- gsub(",", " ", gsName)
    gsName <- gsub(" ", "_", gsName)
    #dir.create(file.path(heatmapDir, gsName), showWarnings = FALSE)
    print(names(gsList)[j])
    print(gsName)
    print(length(gsList[[j]]))
    
    mygenes <- gsList[[j]]
    mygenes <- intersect(mygenes, rownames(expr))
    print(length(mygenes))
    
    # top 50
    overall_heatmap(expr.scaled,
                    annot = ann.sample,
                    degList = degList,
                    mygenes,
                    pv = pvCutoff,
                    fc = fcCutoff,
                    sortIdx = NULL,
                    nbMAX = 50,
                    onlyDEG = FALSE,
                    outFile = file.path(heatmapDir, paste0(gsName, "_overall_heatmap.pdf")))
    
    # ALL  DEG
    overall_heatmap(expr.scaled,
                    annot = ann.sample,
                    degList = degList,
                    mygenes,
                    pv = pvCutoff,
                    fc = fcCutoff,
                    sortIdx = NULL,
                    nbMAX = 1000,
                    onlyDEG = TRUE,
                    outFile = file.path(heatmapDir, paste0(gsName, "_overall_heatmap_DEG.pdf")))
  })
  
})










