# functions
mypeer <- function(expre, HidenFactor = 10) {
  if(!"peer" %in% loadedNamespaces()) library(peer)
  model = PEER()
  PEER_setPhenoMean(model, as.matrix(expre))
  PEER_setNk(model, HidenFactor)
  PEER_getNk(model)
  PEER_setAdd_mean(model, TRUE)
  PEER_setNmax_iterations(model, 10000000)
  PEER_update(model)
  
  factors = PEER_getX(model)
  dim(factors)
  weights = PEER_getW(model)
  dim(weights)
  precision = PEER_getAlpha(model)
  dim(precision)
  residuals = PEER_getResiduals(model)
  dim(residuals)
  #plot(precision)
  # expre_adjust <- mean(rowMeans(expre)) + residuals
  expre_adjust <- residuals
  colnames(expre_adjust) <- colnames(expre)
  rownames(expre_adjust) <- rownames(expre)
  
  return(t(expre_adjust))
}


myQTL <- function(snps, gene, snpspos, genepos) {
  Matrix_eQTL_main(snps = snps,
                   gene = gene,
                   output_file_name = NULL,
                   pvOutputThreshold = 0,
                   useModel = modelLINEAR,
                   errorCovariance = numeric(),
                   verbose = TRUE,
                   output_file_name.cis = NULL,
                   pvOutputThreshold.cis = 1,
                   snpspos = snpspos,
                   genepos = genepos,
                   cisDist = 1e6,
                   pvalue.hist = "qqplot",
                   min.pv.by.genesnp = FALSE,
                   noFDRsaveMemory = FALSE)}



## Load Gene Expression -----------------------------------
library(data.table)
library(peer)
library(MatrixEQTL)

SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt")

# PolyA

PSI_PolyA <- readRDS("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/PolyA_intron-centric_PSI.Rds")
PSI_PolyA[is.na(PSI_PolyA)] <- 0
PSI_PolyA <- PSI_PolyA[, SampInfo_PolyA$donor_id]
stopifnot(identical(colnames(PSI_PolyA), SampInfo_PolyA$donor_id))

# Total

PSI_Total <- readRDS("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/Total_intron-centric_PSI.Rds")
PSI_Total[is.na(PSI_Total)] <- 0
PSI_Total <- PSI_Total[, SampInfo_Total$donor_id]
stopifnot(identical(colnames(PSI_Total), SampInfo_Total$donor_id))

## PEER

PSI_PolyA <- mypeer(expre = t(PSI_PolyA))
PSI_Total <- mypeer(expre = t(PSI_Total))


## load Genotype ------------------------------------------

load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Genotype/snps.RData")
stopifnot(identical(SampInfo_PolyA$donor_id, MatrixEQTL::colnames(snps)))


## Load SNP position --------------------------------------
snpspos <- fread("/mnt/raid61/Personal_data/tangchao/QTL/document/Genotype/SNP_for_QTL/snp_position.txt", header = T, select = 2:4)
setDF(snpspos)


## Load gene position -------------------------------------

library(GenomicRanges)

gr <- as(gsub("\\|s[se]", "", union(row.names(PSI_PolyA), row.names(PSI_Total))), "GRanges")
gr$ID <- union(row.names(PSI_PolyA), row.names(PSI_Total))

genepos <- data.frame(ID = gr$ID, chr = as.character(gsub("chr", "", seqnames(gr))), start = start(gr), end = end(gr))
genepos$chr <- as.character(genepos$chr)

genepos <- subset.data.frame(genepos, chr %in% snpspos$chr)
setDT(genepos)
setkey(genepos, chr, start, end)

PSI_Total <- PSI_Total[row.names(PSI_Total) %in% genepos$ID, ]
PSI_PolyA <- PSI_PolyA[row.names(PSI_PolyA) %in% genepos$ID, ]
setDF(genepos)


## MatrixEQTL ---------------------------------------------

gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(PSI_Total))
PSI_Total_QTL <- myQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos)

save(PSI_Total_QTL, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/analysis/MatrixEQTL/PSI_Total_QTL.RData")

gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(PSI_PolyA))
PSI_PolyA_QTL <- myQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos)

save(PSI_PolyA_QTL, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/analysis/MatrixEQTL/PSI_PolyA_QTL.RData")
