
library(ggplot2)
library(genefilter)
library(data.table)
library(pheatmap)
library(GO.db)
library(org.Hs.eg.db)
library(openxlsx)
library(viridis)
library(gridExtra)
library(limma)
library(msigdbr)
library(parallel)

source("~/Scripts/work/pipelines/enrichment_plot_helper_fct.R")
source("~/Scripts/work/tools/venn_script.r")
source("~/Scripts/work/tools/iwanthue.r")
source("~/Scripts/work/tools/hyperG_parallel.R")
source("~/Scripts/work/tools/jaccardGraph.R")

# Consensus Path DB
#load("~/Research/database/consensus/Consensus_unique08.RData") # load cons2 annotation list
# Misgdb 
#load("~/Research/database/msigdb_human/v7.0/msigdb_v7.0.human.RData")
#dbList[["Consensus"]] <- cons2

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


# ADD CONSENSUS
load(file.path("~/Research/database/consensus_path_db/release34/consensus34.RData"))
dbList[["Consensus"]] <- cons34


############################################
###                                      ###
###               FUNCTION               ###
###                                      ###
############################################

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

toNum <- function(x) return(as.numeric(levels(x))[x])

getDEG <- function(limma, fc, pv, rsign)
{
  limmaFC <- limma$logFC
  limmaPV <- limma$adj.P.Val
  if(rsign == "UP") return(as.character(limma$entrez)[limmaFC >= fc & limmaPV <= pv])
  else if(rsign == "DOWN") return(as.character(limma$entrez)[limmaFC <= -fc & limmaPV <= pv])
  else if(rsign == "BOTH") return(as.character(limma$entrez)[(limmaFC <= -fc | limmaFC >= fc) & limmaPV <= pv])
}

getGenesDESeq2 <- function(deseq, fcCt, pvCt, dsign, adjusted)
{
  deseqFC <- deseq$log2FoldChange
  ifelse(adjusted, deseqPV <- deseq$padj, deseqPV <- deseq$pvalue)
  deseqPV[is.na(deseqPV)] <- 1
  
  deseqUP <- as.character(deseq$ensembl)[deseqFC > fcCt & deseqPV <= pvCt]
  deseqUP <- deseqUP[!deseqUP=="NA"]
  deseqDOWN <- as.character(deseq$ensembl)[deseqFC < -fcCt & deseqPV <= pvCt]
  deseqDOWN <- deseqDOWN[!deseqDOWN=="NA"]
  if(dsign == "up") return(deseqUP)
  else if(dsign == "down") return(deseqDOWN)
  else(return(c(deseqUP, deseqDOWN)))
}


saveFisher <- function(fisherList, outFile)
{
  wb <- createWorkbook()
  for(i in 1:length(fisherList))
  {
    addDataFrame(x = fisherList[[i]], sheet = createSheet(wb = wb, sheetName = names(fisherList)[i]), row.names = FALSE)
  }
  saveWorkbook(wb, outFile)	
}

saveFisherDetails <- function(fisherList, outFile)
{
  fisherList.nrow <- unlist(lapply(fisherList, nrow))
  fisherList <- fisherList[fisherList.nrow != 0]
  fh_info <- data.frame(Sheet = paste("sheet_", seq(1, length(fisherList)), sep = ""),
                        Group = names(fisherList))
  
  toXLSX <- list("info" = fh_info)
  toXLSX <- c(toXLSX, fisherList)
  names(toXLSX)[2:length(toXLSX)] <- paste0("sheet_", seq(1, length(fisherList)))
  
  write.xlsx(toXLSX, outFile, row.names = FALSE, overwrite = TRUE)	
}

FisherFromVenn <- function(venn, db, org.library, univ = NULL)
{
  if(is.null(univ)) univ <- venn$entrez
  groups <- unique(venn$Group)
  
  vList <- lapply(groups, function(i)
    hyperG(venn$entrez[venn$Group == i], geneSets = db, universe = univ, org.library = org.library, cutoff = 1, mincount = 2)
  )
  names(vList) <- groups
  return(vList)
}

FisherFromVennSymbol <- function(venn, db, org.library, univ = NULL)
{
  if(is.null(univ)) univ <- venn[,2]
  univ <- symbol2entrez(univ)
  groups <- unique(venn$Group)
  
  vList <- lapply(groups, function(i){
    entrez.current <- symbol2entrez(venn[venn$Group == i, 2])
    hyperG(entrez.current, geneSets = db, universe = univ, org.library = org.library, cutoff = 1, mincount = 2)
  })
  
  
  names(vList) <- groups
  return(vList)
}


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



############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# PARAMETERS
limmaDir <- file.path("~/Research/XXX/limma_voom/")

fisherDir <- file.path("~/Research/XXX/gsea/Fisher/")
dir.create(fisherDir, recursive = TRUE, showWarnings = FALSE)


###################
# FISHER FROM LIMMA

setwd(limmaDir)
limmaFiles <- list.files(pattern = "_limma.xlsx")
names(limmaFiles) <- gsub("_limma.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

genes <- limmaList[[1]]$entrez
genes <- unique(genes[!is.na(genes)])

pvCutoff <- 0.05
fcCutoff <- 5
limmaList.UP <- lapply(limmaList, function(i) i$entrez[i$adj.P.Val < pvCutoff & i$logFC > fcCutoff])
limmaList.UP <- lapply(limmaList.UP, function(i) unique(i[!is.na(i)]))
limmaList.DOWN <- lapply(limmaList, function(i) i$entrez[i$adj.P.Val < pvCutoff & i$logFC < -fcCutoff])
limmaList.DOWN <- lapply(limmaList.DOWN, function(i) unique(i[!is.na(i)]))
limmaList.DEG <- lapply(limmaList, function(i) i$entrez[i$adj.P.Val < pvCutoff & abs(i$logFC) > fcCutoff])
limmaList.DEG <- lapply(limmaList.DEG, function(i) unique(i[!is.na(i)]))


# Fisher
setwd(fisherDir)
lapply(1:length(dbList), function(i){
  mydb <- dbList[[i]]
  mydb.name <- names(dbList)[i]
  dir.create(mydb.name, showWarnings = FALSE)
  
  fh.UP <- lapply(limmaList.UP, hyperG, geneSets = mydb, universe = genes, org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
  fh.DOWN <- lapply(limmaList.DOWN, hyperG, geneSets = mydb, universe = genes, org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
  fh.DEG <- lapply(limmaList.DEG, hyperG, geneSets = mydb, universe = genes, org.library = "org.Hs.eg.db", cutoff = 1, mincount = 2)
  
  # Save
  for(i in 1:length(fh.UP))
  {
    fhName <- paste(names(fh.UP)[i], "_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx", sep = "")
    write.xlsx(list(UP = fh.UP[[i]], DOWN = fh.DOWN[[i]], DEG = fh.DEG[[i]]), file.path(mydb.name, fhName),
               row.names = FALSE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'), overwrite = TRUE)	
  }
})



#######################################
# HEATMAPS

comp <- names(limmaList)

setwd(fisherDir)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  dir.create(file.path(mydb.name, "heatmaps"), showWarnings = FALSE)
  
  # load fisher results
  fhFiles <- file.path(mydb.name, paste0(comp, "_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx"))
  
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx"), "", basename(fhFiles))
  fhList <- do.call(c, fhList)
  fhList <- fhList[-grep("DEG", names(fhList))]
  
  heatmapPipe(fhList, file.path(mydb.name, "heatmaps", paste0(mydb.name, "_Fisher_heatmap_PV.pdf")),
              nb = 5, TermIdx = 1, pvIdx = 5, qvIdx = NULL, colClust = FALSE)
  heatmapPipe(fhList, file.path(mydb.name, "heatmaps", paste0(mydb.name, "_Fisher_heatmap_QV.pdf")),
              nb = 5, TermIdx = 1, pvIdx = 5, qvIdx = 6, colClust = FALSE)
})

#######################################
# SINGLE BARPLOT

# Load all databases for a given comparison
setwd(fisherDir)
nb <- 50

lapply(comp, function(i){
  
  lapply(1:length(dbList), function(k){
    dir.create(file.path(names(dbList)[k], "barplots"), showWarnings = FALSE)
    
    dbList.sub <- dbList[k]
    fhFiles <- file.path(names(dbList)[k], paste0(i, "_pv", pvCutoff, "_fc", fcCutoff, "_", names(dbList)[k], "_Fisher.xlsx"))
    names(fhFiles) <- paste0(i, "_", names(dbList)[k])
    overallBarplot(fhFiles, nb,
                   outFile = file.path(names(dbList)[k], "barplots", paste0(i, "_", names(dbList)[k], "_top", nb, "_barplot_PV.pdf")),
                   nbIdx = 2, pvIdx = 5, qvIdx = NULL)
    overallBarplot(fhFiles, nb,
                   outFile = file.path(names(dbList)[k], "barplots", paste0(i, "_", names(dbList)[k], "_top", nb, "_barplot_QV.pdf")),
                   nbIdx = 2, pvIdx = 5, qvIdx = 6)
  })
})


#######################################
# OVERALL BARPLOT SELECTED KEYWORDS

# Load all databases for a given comparison
setwd(fisherDir)
dir.create("keywords", showWarnings = FALSE)
dir.create("keywords/barplots", showWarnings = FALSE)
dir.create("keywords/heatmaps", showWarnings = FALSE)
nb <- 20

keywords <- c("metabol", "cytokine", "interleukin")
lapply(keywords, function(kw){
  lapply(comp, function(i){
    dbList.sub <- dbList
    fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_pv", pvCutoff, "_fc", fcCutoff, "_", j, "_Fisher.xlsx")))
    names(fhFiles) <- paste0(i, "_", names(fhFiles))
    overallBarplot(fhFiles, nb,
                   outFile = file.path("keywords/barplots", paste0(i, "_", kw, "_top", nb, "_barplot_PV.pdf")),
                   nbIdx = 2, pvIdx = 5, qvIdx = NULL, kw = kw)
    overallBarplot(fhFiles, nb,
                   outFile = file.path("keywords/barplots", paste0(i, "_", kw, "_top", nb, "_barplot_QV.pdf")),
                   nbIdx = 2, pvIdx = 5, qvIdx = 6, qvreorder = TRUE, kw = kw)
  })
})


#######################################
# OVERALL HEATMAP SELECTED KEYWORDS

setwd(fisherDir)
fhListList <- lapply(comp, function(i){
  dbList.sub <- dbList
  fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_pv", pvCutoff, "_fc", fcCutoff, "_", j, "_Fisher.xlsx")))
  fhList <- lapply(fhFiles, my.read_xlsx)
  fhList <- do.call(c, fhList)
  fhList.up <- do.call(rbind, fhList[grepl("UP$", names(fhList))])
  fhList.down <- do.call(rbind, fhList[grepl("DOWN$", names(fhList))])
  fh <- list(UP = fhList.up, DOWN = fhList.down)
  return(fh)
})
names(fhListList) <- comp
fhList <- do.call(c, fhListList)

lapply(keywords, function(kw){
  heatmapPipe(fh_list = fhList, outFile = file.path("keywords/heatmaps", paste(kw, "_fgsea_heatmap_PV.pdf")),
              nb = 10, TermIdx = 1, pvIdx = 5, kw = kw, gs = NULL, limit = 5, colClust = FALSE)
  heatmapPipe(fh_list = fhList, outFile = file.path("keywords/heatmaps", paste(kw, "_fgsea_heatmap_QV.pdf")),
              nb = 10, TermIdx = 1, pvIdx = 5, qvIdx = 6, kw = kw, gs = NULL, limit = 5, colClust = FALSE)
})


if(FALSE){
  
  
  #######################################
  # OVERALL BARPLOT
  
  # Load all databases for a given comparison
  setwd(fisherDir)
  nb <- 50
  
  lapply(comp, function(i){
    dbList.sub <- dbList
    fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_pv", pvCutoff, "_fc", fcCutoff, "_", j, "_Fisher.xlsx")))
    names(fhFiles) <- paste0(i, "_", names(fhFiles))
    overallBarplot(fhFiles, nb, outName = paste0(i, "_ALL"))
    
    dbList.sub <- dbList[grepl("^C2", names(dbList))]
    fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_pv", pvCutoff, "_fc", fcCutoff, "_", j, "_Fisher.xlsx")))
    names(fhFiles) <- paste0(i, "_", names(fhFiles))
    overallBarplot(fhFiles, nb, outName = paste0(i, "_pathway"))
    
    dbList.sub <- dbList[grepl("^C5", names(dbList))]
    fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_pv", pvCutoff, "_fc", fcCutoff, "_", j, "_Fisher.xlsx")))
    names(fhFiles) <- paste0(i, "_", names(fhFiles))
    overallBarplot(fhFiles, nb, outName = paste0(i, "_GO"))
  })
  
  
  
  #######################################
  # JACCARD GRAPHS
  
  setwd(fisherDir)
  lapply(1:length(dbList), function(i){
    mydb.name <- names(dbList)[i]
    
    fhFiles <- list.files(mydb.name, pattern = paste0("_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx"), full.names = TRUE)
    names(fhFiles) <- gsub(paste0("_pv", pvCutoff, "_fc", fcCutoff, "_", mydb.name, "_Fisher.xlsx"), "", basename(fhFiles))
    fhList.UP <- lapply(fhFiles, read.xlsx, sheet = "UP")
    fhList.DOWN <- lapply(fhFiles, read.xlsx, sheet = "DOWN")
    
    # REMOVE EMPTY
    fhList.UP <- fhList.UP[sapply(fhList.UP, nrow) > 1]
    fhList.DOWN <- fhList.DOWN[sapply(fhList.DOWN, nrow) > 1]
    
    # GET TOP GENE-SETS
    nb.top <- 25
    gsList.UP <- lapply(fhList.UP, getGsFromFisher, top = nb.top)
    gsList.DOWN <- lapply(fhList.DOWN, getGsFromFisher, top = nb.top)
    
    # PLOT JACCARD GRAPH
    dir.create(file.path(names(dbList)[i], "jaccard"), showWarnings = FALSE)
    th <- 0.1
    if(length(gsList.UP)!=0){
      lapply(1:length(gsList.UP), function(j){
        plotJaccard(gsList.UP[[j]],
                    file.path(names(dbList)[i], "jaccard", paste0(names(gsList.UP)[j], "_UP_Jaccard_", th, "_graph.pdf")),
                    threshold = th, score = NULL)
        jaccardHeatmap(gsList.UP[[j]], zeroDiag = TRUE, 
                       file.path(names(dbList)[i], "jaccard", paste0(names(gsList.UP)[j], "_UP_Jaccard_heatmap.pdf")))
      })
    }
    if(length(gsList.DOWN)!=0){
      lapply(1:length(gsList.DOWN), function(j){
        plotJaccard(gsList.DOWN[[j]],
                    file.path(names(dbList)[i], "jaccard", paste0(names(gsList.DOWN)[j], "_DOWN_Jaccard_", th, "_graph.pdf")),
                    threshold = th, score = NULL)
        jaccardHeatmap(gsList.DOWN[[j]], zeroDiag = TRUE, 
                       file.path(names(dbList)[i], "jaccard", paste0(names(gsList.DOWN)[j], "_DOWN_Jaccard_heatmap.pdf")))
      })
    }
  })
  
  
  ########################################
  # FISHER FROM VENN (QUALITATIVE CORRELATION)
  
  vennDir <- file.path("~/Research/XXX/")
  vennFisherDir <- file.path(vennDir, "Fisher")
  dir.create(vennFisherDir, recursive = TRUE, showWarnings = FALSE)
  
  # Load Venn outputs
  setwd(vennDir)
  #vennFiles <- list.files(pattern = "_Venn.xlsx")
  vennFiles <- c("Methy_vs_RNA_vs_Proteome_HM-LM.xlsx")
  names(vennFiles) <- gsub(".xlsx", "", vennFiles)	
  vennList <- lapply(vennFiles, read.xlsx, sheet = "details")	
  
  setwd(vennFisherDir)
  lapply(1:length(dbList), function(k){
    mydb <- dbList[[k]]
    mydb.name <- names(dbList)[k]
    fhList <- lapply(vennList, FisherFromVennSymbol, db = mydb, org.library = "org.Hs.eg.db", univ = NULL)
    
    lapply(1:length(fhList), function(i) 
      saveFisherDetails(fhList[[i]], paste0(names(fhList)[i], "_", mydb.name, "_Fisher.xlsx"))
    )	
  })
}















