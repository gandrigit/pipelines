library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(openxlsx)
library(pheatmap)
library(GO.db)
library(fgsea)
library(viridis)
library(gridExtra)
library(msigdbr)
library(ggplot2)
library(limma)
library(tidyr)
library(clusterProfiler)
library(colorspace)
library(enrichplot)
library(gridExtra)
library(parallel)
library(grid)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

NBCPU <- 8

DCOLORPAL <- "RdBu"
#DCOLORPAL <- "Tropic"
SCOLORPAL <- "Plasma"


#source("~/Scripts/work/tools/fgsea_fct.r")
#source("~/Scripts/work/tools/hyperG.R")
#source("~/Scripts/work/tools/tools.r")
#source("~/Scripts/work/tools/jaccardGraph.R")
#source("~/Scripts/work/tools/map2color.r")
#source("~/Scripts/work/tools/iwanthue.r")
source("~/Scripts/work/pipelines/enrichment_plot_helper_fct.R")


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
#load(file.path("~/Research/database/consensus_path_db/release34/consensus34.RData"))
#dbList[["Consensus"]] <- cons34


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

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
}

entrez2symbol <- function(entrez)
{
  symbol <- mget(as.character(entrez), org.Hs.egSYMBOL, ifnotfound=NA)
  symbol <- unique(unlist(lapply(symbol, function(i) return(i[1]))))
  return(symbol)
}

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

refseq2entrez <- function(refseq){
  entrez <- mget(as.character(refseq), org.Hs.egREFSEQ2EG, ifnotfound=NA)
  entrez <- unique(unlist(lapply(entrez, function(i) return(i[1]))))
  return(entrez)
}

mygsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")

getRankedGenes <- function(degMat)
{
  degMat <- degMat[as.character(degMat$entrez) != "NA", ]
  degMat <- degMat[!duplicated(degMat$entrez), ]
  entrez <- as.character(degMat$entrez) 
  pv <- toNum(degMat$pv)
  fc <- toNum(degMat$fc)
  
  score <- -log10(pv) * fc
  names(score) <- entrez
  
  return(sort(score, decreasing = TRUE))
}

############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################

# SET GENE-SETS LIMITS
dbList <- lapply(dbList, function(db){
  db <- db[sapply(db, length) >= 10]
  return(db)
})

if(FALSE){
  # ENTREZ 2 SYMBOL 
  dbList <- mclapply(dbList, function(db){
    db <- lapply(db, entrez2symbol)
    db <- lapply(db, unique)
    return(db)
  }, mc.cores = NBCPU)
}

# LIST TO TIBBLE
dbTibbles <- lapply(dbList, function(db){
  df <- lapply(1:length(db), function(i){
    return(data.frame(TERM = names(db)[i],
                      geneID = as.character(db[[i]])))
  })
  df <- do.call(rbind, df)
  return(tibble(df))
})


############
# PARAMETERS
limmaDir <- file.path("~/Research/XXX/limma/")

fgseaDir <- file.path("~/Research/XXX/gsea/clusterProfiler_test/")
dir.create(fgseaDir, recursive = TRUE, showWarnings = FALSE)

##################################
# Rank genes based on limma output
setwd(limmaDir)
limmaFiles <- list.files(pattern = "_limma.xlsx", recursive = TRUE)
names(limmaFiles) <- gsub("_limma.xlsx", "", basename(limmaFiles))

limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

rankedGenesList <- lapply(limmaList, getRankedGenes)

if(FALSE){
  ##################################
  # based on TPM (single sample)
  
  countDir <- file.path("~/Research/XXX/RNA/")
  
  setwd(countDir)
  expr <- read.delim("TCGA_RNA_log2TPM.biolinks.txt")
  
  rankedGenesList <- lapply(1:ncol(expr), function(i){
    values <- as.numeric(expr[,i])
    names(values) <- rownames(expr)
    return(sort(values))
  })
}

# Define color
mycolor <- c("#845EC2", "#B0A8B9", "#C34A36", "#008D83")
if(length(mycolor) < length(rankedGenesList)){
  nb.cols <- length(rankedGenesList)
  mycolor <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
}else{
  mycolor <- mycolor[1:length(rankedGenesList)]
}

#########################
# Perform fgsea analysis

setwd(fgseaDir)
lapply(1:length(dbTibbles), function(i){
  mydb <- dbTibbles[[i]]
  mydb.name <- names(dbList)[i]
  
  dir.create(file.path(fgseaDir, mydb.name), showWarnings = FALSE)
  setwd(file.path(fgseaDir, mydb.name))   
  
  gseaList <- lapply(rankedGenesList, GSEA, TERM2GENE = mydb, pvalueCutoff = 1, by = "fgsea",
                     nPermSimple = 10000, maxGSSize = 1000)
  
  #p <- cnetplot(gseaList[[1]], foldChange = rankedGenesList[[1]], colorEdge = TRUE)
  
  # Add symbols
  gseaList.df <- lapply(gseaList, function(i) as.data.frame(i@result))
  
  gseaList.df <- lapply(gseaList.df, function(k){
    k.entrez <- strsplit(k$core_enrichment, split = "\\/")
    k$NB <- sapply(k.entrez, length)
    k.symbol <- lapply(k.entrez, entrez2symbol)
    k.symbol <- lapply(k.symbol, function(j) unique(j[!is.na(j)]))
    k.symbol <- lapply(k.symbol, toString)
    k$leadingEdge.symbol <- unlist(k.symbol)
    k[, colnames(k) != "Description"]
  })
  
  gseaList.df <- lapply(gseaList.df, function(j){
    j$NES[is.na(j$NES)] <- 0
    return(j)
  })
  
  # Divide UP and DOWN
  gseaList.df.UP <- lapply(gseaList.df, function(j) j[j$NES > 0, ])
  names(gseaList.df.UP) <- paste0(names(gseaList.df), ".UP")
  gseaList.df.DOWN <- lapply(gseaList.df, function(j) j[j$NES < 0, ])
  names(gseaList.df.DOWN) <- paste0(names(gseaList.df), ".DOWN")
  gseaList.df <- c(gseaList.df.UP, gseaList.df.DOWN)
  
  # Save
  lapply(1:length(rankedGenesList), function(j)
    write.xlsx(list(UP =  gseaList.df.UP[[j]], DOWN = gseaList.df.DOWN[[j]]),
               paste(names(rankedGenesList)[j], "_", mydb.name, "_fgsea.xlsx", sep = ""),
               rowNames = FALSE, firstRow = T, headerStyle = createStyle(textDecoration = 'bold'), overwrite = TRUE)
  )	
  
  if(FALSE){
    # PLOT
    dir.create(file.path(fgseaDir, mydb.name, "plot"))
    
    nb <- 20
    mygeneSetID <- lapply(gseaList, function(x){
      gs <- x@result$ID
      if(length(gs) > nb) gs <- gs[1:nb]
      return(gs)
    })
    mygeneSetID <- unique(unlist(mygeneSetID))
    if(mydb.name == "H") mygeneSetID <- unique(mydb$TERM)
    
    lapply(mygeneSetID, function(j){
      if(length(gseaList) != 1 & length(gseaList) <= 6){
        p <- mygseaplot_multi(gseaList, geneSetID = j, pvalue_table = TRUE,
                              color = mycolor, ES_geom = "line")
        ggsave(plot = p,
               filename = file.path(fgseaDir, mydb.name, "plot", paste0(j, "_multiplot.pdf")),
               width = 8, height = 6)
      }else{
        dir.create(file.path(fgseaDir, mydb.name, "plot", j))
        lapply(1:length(gseaList), function(x){
          
          p <- mygseaplot2(gseaList[[x]], geneSetID = j, pvalue_table = TRUE,
                           color = "green4", ES_geom = "line")
          ggsave(plot = p,
                 filename = file.path(fgseaDir, mydb.name, "plot", j, paste0(j, "_", names(gseaList)[x], "_plot.pdf")),
                 width = 8, height = 6)
        })
      }
    })
    
  }

})


#######################################
# HEATMAPS

comp <- names(rankedGenesList)

setwd(fgseaDir)
lapply(1:length(dbList), function(i){
  mydb.name <- names(dbList)[i]
  dir.create(file.path(mydb.name, "heatmaps"), showWarnings = FALSE)
  
  # load fisher results
  fhFiles <- file.path(mydb.name, paste0(comp, "_", mydb.name, "_fgsea.xlsx"))
  
  fhList <- lapply(fhFiles, my.read_xlsx)
  names(fhList) <- gsub(paste0("_", mydb.name, "_fgsea.xlsx"), "", fhFiles)
  names(fhList) <- gsub(paste0(mydb.name, "\\/"), "", names(fhList))
  fhList <- do.call(c, fhList)
  
  heatmapPipe(fhList, file.path(mydb.name, "heatmaps", paste0(mydb.name, "_fgsea_heatmap_PV.pdf")),
              nb = 5, TermIdx = 1, pvIdx = 5, qvIdx = NULL, colClust = FALSE)
  heatmapPipe(fhList, file.path(mydb.name, "heatmaps", paste0(mydb.name, "_fgsea_heatmap_QV.pdf")),
              nb = 5, TermIdx = 1, pvIdx = 5, qvIdx = 6, colClust = FALSE)
})

#######################################
# SINGLE BARPLOT

# Load all databases for a given comparison
setwd(fgseaDir)
nb <- 50

lapply(comp, function(i){
  
  lapply(1:length(dbList), function(k){
    dir.create(file.path(names(dbList)[k], "barplots"), showWarnings = FALSE)

    dbList.sub <- dbList[k]
    fhFiles <- file.path(names(dbList)[k], paste0(i, "_", names(dbList)[k], "_fgsea.xlsx"))
    names(fhFiles) <- paste0(i, "_", names(dbList)[k])
    overallBarplot(fhFiles, nb,
                   outFile = file.path(names(dbList)[k], "barplots", paste0(i, "_", names(dbList)[k], "_top", nb, "_barplot_PV.pdf")),
                   nbIdx = 11, pvIdx = 5, qvIdx = NULL)
    overallBarplot(fhFiles, nb,
                   outFile = file.path(names(dbList)[k], "barplots", paste0(i, "_", names(dbList)[k], "_top", nb, "_barplot_QV.pdf")),
                   nbIdx = 11, pvIdx = 5, qvIdx = 6)
  })
})

#######################################
# OVERALL BARPLOT SELECTED KEYWORDS

# Load all databases for a given comparison
setwd(fgseaDir)
dir.create("keywords", showWarnings = FALSE)
dir.create("keywords/barplots", showWarnings = FALSE)
dir.create("keywords/heatmaps", showWarnings = FALSE)
nb <- 20

keywords <- c("metabol", "cytokine", "interleukin")
lapply(keywords, function(kw){
  lapply(comp, function(i){
    dbList.sub <- dbList
    fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_", j, "_fgsea.xlsx")))
    names(fhFiles) <- paste0(i, "_", names(fhFiles))
    overallBarplot(fhFiles, nb,
                   outFile = file.path("keywords/barplots", paste0(i, "_", kw, "_top", nb, "_barplot_PV.pdf")),
                   nbIdx = 11, pvIdx = 5, qvIdx = NULL, kw = kw)
    overallBarplot(fhFiles, nb,
                   outFile = file.path("keywords/barplots", paste0(i, "_", kw, "_top", nb, "_barplot_QV.pdf")),
                   nbIdx = 11, pvIdx = 5, qvIdx = 6, qvreorder = TRUE, kw = kw)
  })
})

#######################################
# OVERALL HEATMAP SELECTED KEYWORDS

setwd(fgseaDir)
fhListList <- lapply(comp, function(i){
  dbList.sub <- dbList
  fhFiles <- sapply(names(dbList.sub), function(j) file.path(j, paste0(i, "_", j, "_fgsea.xlsx")))
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





