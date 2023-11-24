
library(limma)
library(edgeR)
library(openxlsx)
library(UpSetR)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)

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

ensembl2entrez <- function(ensembl)
{
  entrez <- mget(as.character(ensembl), org.Hs.egENSEMBL2EG, ifnotfound=NA)
  entrez <- lapply(entrez, function(i) return(i[1]))
  return(unlist(entrez))
}

symbol2entrez <- function(symbol)
{
  entrez <- mget(as.character(symbol), org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrez <- unlist(lapply(entrez, function(i) return(i[1])))
  entrez <- unique(entrez[!is.na(entrez)])
  return(entrez)
}

entrez2name <- function(entrez)
{
  gn <- mget(as.character(entrez), org.Hs.egGENENAME, ifnotfound=NA)
  gn <- unlist(lapply(gn, function(i) return(i[1])))
  return(gn)
}


getDEG <- function(limma, pv, fc, dsign = "both", adjusted = TRUE)
{
  limmaFC <- limma$logFC
  if(class(limmaFC)=="factor") limmaFC <- as.numeric(levels(limmaFC))[limmaFC]
  ifelse(adjusted, limmaPV <- limma$adj.P.Val, limmaPV <- limma$P.Value)
  if(class(limmaPV)=="factor") limmaPV <- as.numeric(levels(limmaPV))[limmaPV]
  
  if(dsign == "up") return(as.character(limma$entrez)[limmaFC > fc & limmaPV < pv])
  else if(dsign == "down") return(as.character(limma$entrez)[limmaFC < -fc & limmaPV < pv])
  else if(dsign == "both") return(as.character(limma$entrez)[abs(limmaFC) > fc & limmaPV < pv])
}

getNbDEG <- function(limma, pv, fc, adjusted = FALSE)
{
  limmaFC <- limma$logFC
  if(class(limmaFC)=="factor") limmaFC <- as.numeric(levels(limmaFC))[limmaFC]
  ifelse(adjusted, limmaPV <- limma$adj.P.Val, limmaPV <- limma$P.Value)
  if(class(limmaPV)=="factor") limmaPV <- as.numeric(levels(limmaPV))[limmaPV]
  
  return(sum(abs(limmaFC)>= fc & limmaPV <= pv, na.rm = TRUE))
}


ggVolcano <- function(limma, pvalue = 0, genes = NULL, outFile){
  ggmat <- data.frame(X = limma$logFC, Y = -log10(limma$adj.P.Val), GENE = limma$symbol)
  ggmat$isSignif <- "no"
  ggmat$isSignif[limma$adj.P.Val < pvalue] <- "yes"
  ggmat$toShow <- "no"
  if(!is.null(genes)){
    ggmat$toShow <- "no"
    ggmat$toShow[match(genes, ggmat$GENE)] <- "yes"
    
    ggmat$isSignif <- "no"
    ggmat$isSignif[match(genes, ggmat$GENE)] <- "yes"
  }
  
  p <- ggplot(ggmat, aes(x=X, y=Y, color = isSignif, label = GENE))
  p <- p + geom_point(size = 2, alpha = 0.5, data = subset(ggmat, isSignif == "no"))
  p <- p + geom_point(size = 2, alpha = 0.5, data = subset(ggmat, isSignif == "yes"))
  p <- p + theme_bw(base_size = 20)
  p <- p + theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18))
  p <- p + xlab("log2 Fold Change")+ ylab("-log10 Pvalue")
  p <- p + theme(legend.position="none")
  p <- p + scale_color_manual(values=c(yes = "red2", no = "grey"))
  p <- p + geom_text_repel(data = subset(ggmat, toShow == "yes"), color = "black", force = 2)
  
  pdf(outFile)
  plot(p)
  dev.off()
}

ggVolcanoPV <- function(limma, pvalue = 0, genes = NULL, outFile){
  ggmat <- data.frame(X = limma$logFC, Y = -log10(limma$P.Value), GENE = limma$symbol, AVG = limma$AveExpr)
  ggmat$isSignif <- "no"
  ggmat$isSignif[limma$P.Value < pvalue] <- "yes"
  ggmat$toShow <- "no"
  if(!is.null(genes)){
    ggmat$toShow <- "no"
    ggmat$toShow[match(genes, ggmat$GENE)] <- "yes"
    
    ggmat$isSignif <- "no"
    ggmat$isSignif[match(genes, ggmat$GENE)] <- "yes"
  }
  
  p <- ggplot(ggmat, aes(x=X, y=Y, color = AVG, label = GENE))
  p <- p + geom_point(size = 2, alpha = 0.5) + scale_color_viridis_c()
  p <- p + geom_point(size = 2, alpha = 0.5, shape = 21, color = "black", data = subset(ggmat, toShow == "yes"))
  if(pvalue != 0) p <- p + geom_hline(yintercept=-log10(pvalue), linetype="dashed", color = "black")
  p <- p + theme_bw(base_size = 20)
  p <- p + theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18))
  p <- p + xlab("log2 Fold Change")+ ylab("-log10 Pvalue")
  #p <- p + theme(legend.position="none")
  p <- p + geom_text_repel(data = subset(ggmat, toShow == "yes"), color = "black", box.padding = 0.5, max.overlaps = Inf)
  
  pdf(outFile)
  plot(p)
  dev.off()
}

ggVolcanoQV <- function(limma, pvalue = 0, genes = NULL, outFile){
  ggmat <- data.frame(X = limma$logFC, Y = -log10(limma$adj.P.Val), GENE = limma$symbol, AVG = limma$AveExpr)
  ggmat$isSignif <- "no"
  ggmat$isSignif[limma$adj.P.Val < pvalue] <- "yes"
  ggmat$toShow <- "no"
  if(!is.null(genes)){
    ggmat$toShow <- "no"
    ggmat$toShow[match(genes, ggmat$GENE)] <- "yes"
    
    ggmat$isSignif <- "no"
    ggmat$isSignif[match(genes, ggmat$GENE)] <- "yes"
  }
  
  p <- ggplot(ggmat, aes(x=X, y=Y, color = AVG, label = GENE))
  p <- p + geom_point(size = 2, alpha = 0.5) + scale_color_viridis_c()
  p <- p + geom_point(size = 2, alpha = 0.5, shape = 21, color = "black", data = subset(ggmat, toShow == "yes"))
  if(pvalue != 0) p <- p + geom_hline(yintercept=-log10(pvalue), linetype="dashed", color = "black")
  p <- p + theme_bw(base_size = 20)
  p <- p + theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18))
  p <- p + xlab("log2 Fold Change")+ ylab("-log10 Pvalue")
  #p <- p + theme(legend.position="none")
  p <- p + geom_text_repel(data = subset(ggmat, toShow == "yes"), color = "black", box.padding = 0.5, max.overlaps = Inf)
  
  pdf(outFile)
  plot(p)
  dev.off()
}

ggScatterplot <- function(limmaX, limmaY, top = NULL, genes = NULL, outFile){
  common <- intersect(limmaX$entrez, limmaY$entrez)
  limmaX <- limmaX[match(common, limmaX$entrez), ]
  limmaY <- limmaY[match(common, limmaY$entrez), ]
  
  ggmat <- data.frame(X = limmaX$logFC, Y = limmaY$logFC, GENE = limmaX$symbol)
  ggmat$AVG <- rowMeans(ggmat[, c("X", "Y")])
  ggmat <- ggmat[order(ggmat$AVG), ]
  
  ggmat$toShow <- "no"
  if(!is.null(genes)){
    ggmat$toShow <- "no"
    ggmat$toShow[match(genes, ggmat$GENE)] <- "yes"
  }else if(!is.null(top)){
    idx.up <- seq(1, top)
    idx.down <- seq(nrow(ggmat)-top+1, nrow(ggmat))
    ggmat$toShow[c(idx.up, idx.down)] <- "yes"
  }
  
  
  mycor <- cor.test(ggmat$X, ggmat$Y, method = "pearson")
  mytitle <- paste0("R=", round(mycor$estimate, digits = 2), "; pvalue=", signif(mycor$p.value, digits = 3))
  
  p <- ggplot(ggmat, aes(x=X, y=Y, color = AVG, label = GENE))
  p <- p + geom_point(size = 3, alpha = 1)
  p <- p + scale_color_gradient2(midpoint=0,  low="deepskyblue", mid="snow", high="orangered")
  p <- p + geom_point(size = 3, alpha = 1, shape = 21, color = "black", data = subset(ggmat, toShow == "yes"))
  p <- p + theme_bw(base_size = 20)
  p <- p + geom_smooth(method='lm', color = "black")
  p <- p + theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18))
  p <- p + xlab("log2 Fold Change")+ ylab("log2 Fold Change") + ggtitle(mytitle)
  #p <- p + theme(legend.position="none")
  p <- p + geom_text_repel(data = subset(ggmat, toShow == "yes"), color = "black", box.padding = 0.5)
  
  pdf(outFile)
  plot(p)
  dev.off()
}


getNbMat <- function(limma, pvs, fcs, sign="ALL"){
  nbMat <- c()
  for(i in pvs){
    for(j in fcs){
      if(sign == "ALL"){
        nb <- sum(limma$adj.P.Val < i & abs(limma$logFC) > j)
      }else if(sign == "UP"){
        nb <- sum(limma$adj.P.Val < i & limma$logFC > j)
      }else if(sign == "DOWN"){
        nb <- sum(limma$adj.P.Val < i & limma$logFC < -j)
      }else(return(NA))
      nbMat <- rbind(nbMat, c(i, j, nb))
    }
  }
  colnames(nbMat) <- c("Pvalue", "logFC", "NB")
  
  nbMat <- as.data.frame(nbMat)
  nbMat$Pvalue <- as.character(nbMat$Pvalue)
  nbMat$Sign <- sign
  return(nbMat)
}

ggDEG <- function(limmaList, pvCutoff = 0.05, fcCutoff = 0, adjusted = TRUE, outFile){
  
  if(adjusted){
    nbUP <- sapply(limmaList, function(i) sum(i$adj.P.Val < pvCutoff & i$logFC > fcCutoff))
    nbDOWN <- sapply(limmaList, function(i) sum(i$adj.P.Val < pvCutoff & i$logFC < -fcCutoff))
  }else{
    nbUP <- sapply(limmaList, function(i) sum(i$P.Value < pvCutoff & i$logFC > fcCutoff))
    nbDOWN <- sapply(limmaList, function(i) sum(i$P.Value < pvCutoff & i$logFC < -fcCutoff))
  }
  nbTOT <- nbUP + nbDOWN
  
  ggmat <- data.frame(COMP = c(names(nbUP), names(nbDOWN)),
                      SIGN = c(rep("UP", length(nbUP)), rep("DOWN", length(nbDOWN))),
                      NB = c(nbUP, -nbDOWN))
  ggmat$SIGN = factor(ggmat$SIGN, levels = c("UP", "DOWN"))
  ggmat$COMP <- factor(ggmat$COMP, levels = names(nbTOT)[order(nbTOT)])
  
  p <- ggplot(data = ggmat, aes(x = COMP, y = SIGN)) +
    geom_tile(aes(fill = NB)) +
    coord_flip() +
    geom_text(aes(label = abs(NB))) +
    scale_fill_gradient2(low = "skyblue", mid = "snow", high = "tomato") +
    theme_minimal() + 
    theme(legend.position = "none") + 
    xlab("") + ylab("")
  
  myheight <- 2 + 0.25 * length(limmaList)
  ggsave(plot = p, filename = outFile,
         height = myheight, width = 4)
}


log2df <- function(logM)
{
  reads.tot <- logM$V2[logM$V1 == "                          Number of input reads |"]
  mapped <- logM$V2[logM$V1 == "                   Uniquely mapped reads number |"]
  mapped.pc <- logM$V2[logM$V1 == "                        Uniquely mapped reads % |"]
  multi <- logM$V2[logM$V1 == "        Number of reads mapped to multiple loci |"]
  multi.pc <- logM$V2[logM$V1 == "             % of reads mapped to multiple loci |"]
  unmapped.mismatches <- logM$V2[logM$V1 == "       % of reads unmapped: too many mismatches |"]
  unmapped.short <- logM$V2[logM$V1 == "                 % of reads unmapped: too short |"]
  unmapped.other <- logM$V2[logM$V1 == "                     % of reads unmapped: other |"]
  chimeric <- logM$V2[logM$V1 == "                            % of chimeric reads |"]
  
  return(data.frame(Reads.TOTAL = reads.tot,
                    Mapped = mapped,
                    Mapped.PC = mapped.pc,
                    MultipleLoci = multi,
                    MultipleLoci.PC = multi.pc,
                    UnMapped.mismatches.PC = unmapped.mismatches,
                    UnMapped.tooShort.PC = unmapped.short,
                    UnMapped.other.pc = unmapped.other,
                    Chimeric.PC = chimeric)
  )
  
}


############################################
###                                      ###
###                 MAIN                 ###
###                                      ###
############################################


# PARAMETERS
quantDir <- file.path("~/cluster/cluster/XXX/quant/")

docDir <- file.path("~/Research/XXX/doc/")

countDir <- file.path("~/Research/XXX/count/")
dir.create(countDir, recursive = TRUE, showWarnings = FALSE)

limmaDir <- file.path("~/Research/XXX/limma_voom/")
dir.create(limmaDir, recursive = TRUE, showWarnings = FALSE)

##############
# LOAD SAMPLES


if(FALSE){
  #################
  # LOG 2 DATAFRAME
  setwd(quantDir)
  logFiles <- list.files(pattern = "_Log.final.out$", recursive = TRUE)
  names(logFiles) <- dirname(logFiles)
  
  logList <- lapply(logFiles, read.delim, header = FALSE)
  
  logDF <- lapply(logList, log2df)
  logDF <- do.call(rbind, logDF)
  
  setwd(countDir)
  write.xlsx(logDF, "RNA_Log.final.out.xlsx", row.names = TRUE, overwrite = TRUE)
  
  
  #################
  # LOAD COUNT DATA
  setwd(quantDir)
  countFiles <- list.files(pattern = "_ReadsPerGene.out.tab", recursive = TRUE)
  names(countFiles) <- dirname(countFiles)
  countList <- lapply(countFiles, read.delim, skip = 4, header = FALSE, stringsAsFactors = FALSE)
  
  
  # Create geneMat
  ensID <- countList[[1]]$V1
  entrezID <- ensembl2entrez(ensID)
  entrezID[is.na(entrezID)] <- "NA"
  symbol <- entrez2symbol(entrezID)
  geneName <- entrez2name(entrezID)
  
  geneMat <- data.frame(ensembl = ensID, entrez = entrezID, symbol = symbol, gene.name = geneName)
  
  countList <- lapply(countList, function(i) return(i[,2]))
  countMat <- do.call(cbind, countList)
  
  setwd(countDir)
  write.table(countMat, "count.raw.txt", sep = "\t", quote = FALSE)
}

setwd(countDir)
countMat <- read.delim("count.raw.txt")
#rownames(countMat) <- ensID

# Min row CPM >= 2
#logcpmMat <- cpm(countMat, log = TRUE)
#keep.exprs <- rowSums(logcpmMat > 0) >= 2

# count != 0 in at least 2 samples
nonZero <- apply(countMat, 1, function(i) sum(i != 0))
keep.exprs <- nonZero >= 2

# ANNOTATIONS
setwd(docDir)
ann.sample <- read.xlsx("RNA_sample_annotation.xlsx", sheet = 1)
ann.sample <- ann.sample[match(colnames(countMat), ann.sample$SAMPLE), ]

# BUILD DGE
dge <- DGEList(count = countMat, genes = geneMat)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge, method = "TMM") # normalization

m <- dge$counts
rownames(m) <- dge$genes$ensembl


# SAVE COUNT
setwd(countDir)
write.table(m, "count.txt", sep = "\t", quote = FALSE)
#write.table(cpm(m, normalized.lib.sizes = TRUE, log = TRUE), "logCPM.txt", sep = "\t", quote = FALSE)
write.table(log2(cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)+1), "logCPM.txt", sep = "\t", quote = FALSE)



###################################################################################
# LIMMA 

# PLEASE CHECK "groups", "design" and "contrast.matrix" variables

#######################
# Differential analysis
groups <- ann.sample$GROUP

design <- model.matrix(~0+groups) # non-paired
colnames(design) <- gsub("groups", "", colnames(design))
v <- voom(dge, design)

# UNPAIRED
fit <- lmFit(v, design, method = "robust")# robust or not 

fit <- eBayes(fit)
plotSA(fit)

tt <- topTable(fit, adjust="BH", number=nrow(dge), p.value = 1)

contrast.matrix <- makeContrasts(
  X-Y,
  levels=design
)
contrast.name <- gsub(" ", "", colnames(contrast.matrix))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)	

resList <- lapply(1:ncol(contrast.matrix), function(i)
  topTable(fit2, coef=i, adjust="BH", number=nrow(dge), sort.by="P", p.value = 1)
)
names(resList) <- contrast.name

# SAVE results
setwd(limmaDir)
lapply(1:length(resList), function(i)
  write.xlsx(resList[[i]], paste(names(resList)[i], "_limma.xlsx", sep = ""),
             firstRow = T, headerStyle = createStyle(textDecoration = 'bold'), overwrite = TRUE)
)



###################################################################################
# VOLCANO PLOTS


# Load limma outputs
setwd(limmaDir)
limmaFiles <- list.files(pattern = "_limma.xlsx")
names(limmaFiles) <- gsub("_limma.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

limmaList <- lapply(limmaList, function(i) i[!is.na(i$symbol), ])

lapply(1:length(limmaList), function(i){
  ggVolcanoPV(limmaList[[i]], pvalue = 0.05, genes = head(limmaList[[i]]$symbol, 25),
            outFile = paste0(names(limmaList)[i], "_volcano_PV.pdf"))
  ggVolcanoQV(limmaList[[i]], pvalue = 0.05, genes = head(limmaList[[i]]$symbol, 25),
            outFile = paste0(names(limmaList)[i], "_volcano_QV.pdf"))
})



###################################################################################
# NUMBER OF DEG PLOT

setwd(limmaDir)
limmaFiles <- list.files(pattern = "_limma.xlsx")
names(limmaFiles) <- gsub("_limma.xlsx", "", limmaFiles)
limmaList <- lapply(limmaFiles, read.xlsx, sheet = 1)

pvs <- c(0.1, 0.05, 0.01, 0.005, 0.001)
fcs <- seq(0,3, by = 0.5)

for(i in 1:length(limmaList)){
  limmaMat <- limmaList[[i]]
  
  # plot
  ggmat.ALL <- getNbMat(limmaMat, pvs, fcs, sign="ALL")
  ggmat.UP <- getNbMat(limmaMat, pvs, fcs, sign="UP")
  ggmat.DOWN <- getNbMat(limmaMat, pvs, fcs, sign="DOWN")
  
  ggmat <- do.call(rbind,list(ggmat.ALL, ggmat.UP, ggmat.DOWN))
  
  p <- ggplot(data=ggmat, aes(x=logFC, y=NB, group=Pvalue, colour=Pvalue))
  p <- p + theme_bw()
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + facet_grid(. ~ Sign)
  #p <- p + scale_y_log10()
  p <- p + ylab("Number of DEG")
  pdf(paste0(names(limmaList)[i], "_nb_diff.pdf"), width = 10, height = 5)
  plot(p)
  dev.off()
  
}


###################################################################################
# NUMBER OF DEG HEATMAP

ggDEG(limmaList, pvCutoff = 0.05, fcCutoff = 0, adjusted = TRUE,
      outFile = "NB_DEG_QV_heatmap.pdf")

ggDEG(limmaList, pvCutoff = 0.05, fcCutoff = 0, adjusted = FALSE,
      outFile = "NB_DEG_PV_heatmap.pdf")

