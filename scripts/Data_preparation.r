# Data preparation
library(data.table)
## Sample information
SampInfo <- fread("/mnt/raid61/Personal_data/tangchao/QTL/document/SampleInfo/EGAD00001002671_SampInfo.txt")
SampInfo_mRNA <- SampInfo[MOLECULE != "total RNA", ]
SampInfo_tRNA <- SampInfo[donor_id %in% SampInfo_mRNA$donor_id & MOLECULE == "total RNA", ]
fwrite(SampInfo_mRNA, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt", sep = "\t")
fwrite(SampInfo_tRNA, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt", sep = "\t")

## Genotype
snp <- fread("/mnt/raid61/Personal_data/tangchao/IR6/document/genotype/SNP_for_QTL/snp.t.txt")
snp <- data.frame(snp[, -1], row.names = snp[[1]])
snp <- snp[, SampInfo_mRNA$donor_id]

## Outliers in genotype. Minor allele frequency filtering.
maf = rowMeans(snp, na.rm = TRUE)/2
maf = pmin(maf, 1-maf)
snp <- snp[maf > 0.05, ]

fwrite(snp, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Genotype/SNP_Matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE)
system("bgzip /mnt/raid61/Personal_data/tangchao/ScientificData/data/Genotype/SNP_Matrix.txt")

library(MatrixEQTL)
snps = SlicedData$new()
snps$CreateFromMatrix(as.matrix(snp))
save(snps, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Genotype/snps.RData")

## Phenotype
### Gene 
setwd("/mnt/raid62/Chen_Cell_2016/RSEM/EGAD00001002671")
files_PolyA <- paste(SampInfo_mRNA$ID , ".genes.results", sep = "")
stopifnot(all(file.exists(files_PolyA)))
Gene_PolyA <- lapply(files_PolyA, fread)
expected_count_PolyA <- do.call(cbind, lapply(Gene_PolyA, function(x) x[,.(expected_count)]))
TPM_PolyA <- do.call(cbind, lapply(Gene_PolyA, function(x) x[,.(TPM)]))
colnames(expected_count_PolyA) <- colnames(TPM_PolyA) <- gsub(".genes.results", "", basename(files_PolyA))


files_Total <- paste(SampInfo_tRNA$ID , ".genes.results", sep = "")
stopifnot(all(file.exists(files_Total)))
Gene_Total <- lapply(files_Total, fread)
expected_count_Total <- do.call(cbind, lapply(Gene_Total, function(x) x[,.(expected_count)]))
TPM_Total <- do.call(cbind, lapply(Gene_Total, function(x) x[,.(TPM)]))
colnames(expected_count_Total) <- colnames(TPM_Total) <- gsub(".genes.results", "", basename(files_Total))


suppressPackageStartupMessages(library(rtracklayer))
gtf <- rtracklayer::readGFF("/mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/gencode.v30lift37.annotation.gtf")
setDT(gtf)
gtf <- gtf[type == "gene", ]
identical(Gene_Total[[1]][[1]], gtf[, paste(gene_id, gene_name, sep="_")])
setkey(gtf, gene_id)

identical(Gene_PolyA[[1]][[1]], gtf[,paste(gene_id, gene_name, sep="_")])
expected_count_PolyA <- cbind(gtf[, .(seqid, start, end, strand, gene_id, gene_type, gene_name)], expected_count_PolyA)
TPM_PolyA <- cbind(gtf[,.(seqid, start, end, strand, gene_id, gene_type, gene_name)], TPM_PolyA)

identical(Gene_Total[[1]][[1]], gtf[,paste(gene_id, gene_name, sep="_")])
expected_count_Total <- cbind(gtf[, .(seqid, start, end, strand, gene_id, gene_type, gene_name)], expected_count_Total)
TPM_Total <- cbind(gtf[,.(seqid, start, end, strand, gene_id, gene_type, gene_name)], TPM_Total)

save(expected_count_PolyA, TPM_PolyA, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/PolyA_RNA_raw_RSEM.RData")
save(expected_count_Total, TPM_Total, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/Total_RNA_raw_RSEM.RData")


### Transcript

files_PolyA <- paste(SampInfo_mRNA$ID , ".isoforms.results", sep = "")
stopifnot(all(file.exists(files_PolyA)))
Gene_PolyA <- lapply(files_PolyA, fread)
expected_count_PolyA <- do.call(cbind, lapply(Gene_PolyA, function(x) x[,.(expected_count)]))
TPM_PolyA <- do.call(cbind, lapply(Gene_PolyA, function(x) x[,.(TPM)]))
colnames(expected_count_PolyA) <- colnames(TPM_PolyA) <- gsub(".isoforms.results", "", basename(files_PolyA))

files_Total <- paste(SampInfo_tRNA$ID , ".isoforms.results", sep = "")
stopifnot(all(file.exists(files_Total)))
Gene_Total <- lapply(files_Total, fread)
expected_count_Total <- do.call(cbind, lapply(Gene_Total, function(x) x[,.(expected_count)]))
TPM_Total <- do.call(cbind, lapply(Gene_Total, function(x) x[,.(TPM)]))
colnames(expected_count_Total) <- colnames(TPM_Total) <- gsub(".isoforms.results", "", basename(files_Total))


gtf <- rtracklayer::readGFF("/mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/gencode.v30lift37.annotation.gtf")
setDT(gtf)
gtf <- gtf[type == "transcript", ]
gtf[, tx_id:=substring(transcript_id, 1, 15)]
gtf <- gtf[match(substring(Gene_PolyA[[1]][[1]], 1, 15), tx_id), ]

stopifnot(identical(substring(Gene_PolyA[[1]][[1]], 1, 15), gtf[, tx_id]))
expected_count_PolyA <- cbind(gtf[,.(seqid, start, end, strand, transcript_id, transcript_type, transcript_name)], expected_count_PolyA)
TPM_PolyA <- cbind(gtf[,.(seqid, start, end, strand, transcript_id, transcript_type, transcript_name)], TPM_PolyA)

stopifnot(identical(substring(Gene_Total[[1]][[1]], 1, 15), gtf[, tx_id]))
expected_count_Total <- cbind(gtf[,.(seqid, start, end, strand, transcript_id, transcript_type, transcript_name)], expected_count_Total)
TPM_Total <- cbind(gtf[,.(seqid, start, end, strand, transcript_id, transcript_type, transcript_name)], TPM_Total)

save(expected_count_PolyA, TPM_PolyA, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/PolyA_RNA_raw_RSEM.RData")
save(expected_count_Total, TPM_Total, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/Total_RNA_raw_RSEM.RData")


### Junction

library(GenomicAlignments)

setwd("/mnt/raid62/Chen_Cell_2016/EGAD00001002671_bam")
files <- paste(SampInfo_mRNA$ID , ".SJ.out.tab", sep = "")
stopifnot(all(file.exists(files)))

mclapply(files, function(x) {
  sj <- readSTARJunctions(x)[, 3]
  sj <- sj[sj$um_reads > 0 & seqnames(sj) %in% paste0("chr", c(1:22, "X", "Y")), ]
  sj <- as.data.table(sj)
  setnames(sj, "um_reads", gsub(".SJ.out.tab", "", basename(x)))
  setkey(sj, seqnames, start, end, width, strand)
  return(sj)
}, mc.cores = 16) -> SJ_list

all.sj <- unique(do.call(rbind, lapply(SJ_list, function(x) x[, 1:5])))

lapply(SJ_list, function(x) {
  x[all.sj, ]
}) -> SJ_list

SJ_tab <- cbind(all.sj, do.call(cbind, lapply(SJ_list, function(x) x[, 6])))
save(SJ_tab, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/PolyA_RNA_Raw_SJs.RData")


files <- paste(SampInfo_tRNA$ID , ".SJ.out.tab", sep = "")
stopifnot(all(file.exists(files)))

mclapply(files, function(x) {
  sj <- readSTARJunctions(x)[, 3]
  sj <- sj[sj$um_reads > 0 & seqnames(sj) %in% paste0("chr", c(1:22, "X", "Y")), ]
  sj <- as.data.table(sj)
  setnames(sj, "um_reads", gsub(".SJ.out.tab", "", basename(x)))
  setkey(sj, seqnames, start, end, width, strand)
  return(sj)
}, mc.cores = 16) -> SJ_list

all.sj <- unique(do.call(rbind, lapply(SJ_list, function(x) x[, 1:5])))

lapply(SJ_list, function(x) {
  x[all.sj, ]
}) -> SJ_list

SJ_tab <- cbind(all.sj, do.call(cbind, lapply(SJ_list, function(x) x[, 6])))
save(SJ_tab, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/Total_RNA_Raw_SJs.RData")

