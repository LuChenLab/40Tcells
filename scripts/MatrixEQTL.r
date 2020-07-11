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

# RSEM

RSEM_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/PolyA_RNA_TPM.txt.gz")
RSEM_PolyA <- data.frame(RSEM_PolyA[, -1], row.names = RSEM_PolyA[[1]])
stopifnot(identical(colnames(RSEM_PolyA), SampInfo_PolyA$donor_id))

RSEM_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/Total_RNA_TPM.txt.gz")
RSEM_Total <- data.frame(RSEM_Total[, -1], row.names = RSEM_Total[[1]])
stopifnot(identical(colnames(RSEM_Total), SampInfo_Total$donor_id))

# HTSeq

library(GenomicAlignments)
library(DESeq2)
load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/se.RData")

dds_PolyA <- DESeqDataSet(se_PolyA, design = ~ 1)
dds_PolyA <- dds_PolyA[rowSums(counts(dds_PolyA) >= 10) >= .2*ncol(dds_PolyA), ]
dds_PolyA <- estimateSizeFactors(dds_PolyA)
dds_PolyA <- estimateDispersions(dds_PolyA)

dds_Total <- DESeqDataSet(se_Total, design = ~ 1)
dds_Total <- dds_Total[rowSums(counts(dds_Total) >= 10) >= .2*ncol(dds_Total), ]
dds_Total <- estimateSizeFactors(dds_Total)
dds_Total <- estimateDispersions(dds_Total)

HTSeq_PolyA <- counts(dds_PolyA, normalized = TRUE)
HTSeq_Total <- counts(dds_Total, normalized = TRUE)

colnames(HTSeq_PolyA) <- gsub(".Aligned.sortedByCoord.out.bam", "", colnames(HTSeq_PolyA))
colnames(HTSeq_Total) <- gsub(".Aligned.sortedByCoord.out.bam", "", colnames(HTSeq_Total))
if(identical(colnames(HTSeq_PolyA), SampInfo_PolyA$ID)) colnames(HTSeq_PolyA) <- SampInfo_PolyA$donor_id
if(identical(colnames(HTSeq_Total), SampInfo_Total$ID)) colnames(HTSeq_Total) <- SampInfo_Total$donor_id
stopifnot(identical(colnames(HTSeq_PolyA), SampInfo_PolyA$donor_id))
stopifnot(identical(colnames(HTSeq_Total), SampInfo_Total$donor_id))

identical(SampInfo_Total$donor_id, SampInfo_PolyA$donor_id)

## PEER
RSEM_PolyA <- mypeer(expre = t(log2(1 + RSEM_PolyA)))
RSEM_Total <- mypeer(expre = t(log2(1 + RSEM_Total)))

HTSeq_PolyA <- mypeer(expre = t(log2(1 + HTSeq_PolyA)))
HTSeq_Total <- mypeer(expre = t(log2(1 + HTSeq_Total)))


## load Genotype ------------------------------------------

load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Genotype/snps.RData")
stopifnot(identical(SampInfo_PolyA$donor_id, MatrixEQTL::colnames(snps)))

## Load SNP position --------------------------------------
snpspos <- fread("/mnt/raid61/Personal_data/tangchao/QTL/document/Genotype/SNP_for_QTL/snp_position.txt", header = T, select = 2:4)
setDF(snpspos)


## Load gene position -------------------------------------
gtf <- rtracklayer::readGFF("/mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/gencode.v30lift37.annotation.gtf")
setDT(gtf)
gtf <- gtf[type == "gene", .(seqid, start, end, strand, gene_id, gene_type, gene_name)]
gtf <- gtf[,.(gene_id, seqid, start, end)]
gtf[, gene_id:=as.character(gene_id)]
gtf <- gtf[seqid %in% paste("chr", c(1:22, "M"), sep = ""),]
genepos <- as.data.frame(gtf)
colnames(genepos) <- c("ID", "chr", "start", "end")
genepos$chr <- gsub("chr", "", genepos$chr)



## MatrixEQTL ---------------------------------------------

gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(RSEM_PolyA))
RSEM_PolyA_QTL <- myQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos)


gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(RSEM_Total))
RSEM_Total_QTL <- myQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos)


gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(HTSeq_PolyA))
HTSeq_PolyA_QTL <- myQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos)


gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(HTSeq_Total))
HTSeq_Total_QTL <- myQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos)

save(RSEM_PolyA_QTL, RSEM_Total_QTL, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/analysis/MatrixEQTL/RSEM_eQTL.RData")
save(HTSeq_PolyA_QTL, HTSeq_Total_QTL, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/analysis/MatrixEQTL/HTSeq_eQTL.RData")

