# Gene
library(data.table)
## PloyA RNA
SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/PolyA_RNA_raw_RSEM.RData")
Indi1 <- SampInfo_PolyA$ID
Index <- rowSums(expected_count_PolyA[, ..Indi1] >= 6) >= 0.2*length(Indi1) & rowSums(TPM_PolyA[, ..Indi1] >= 0.1) >= 0.2*length(Indi1)
ExpMat <- data.frame(TPM_PolyA[Index, ..Indi1])
row.names(ExpMat) <- TPM_PolyA[Index, gene_id]
stopifnot(identical(colnames(ExpMat), SampInfo_PolyA$ID))
colnames(ExpMat) <- SampInfo_PolyA$donor_id
TPM_PolyA <- ExpMat

## Total RNA
SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt", sep = "\t")
load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/Total_RNA_raw_RSEM.RData")
Indi1 <- SampInfo_Total$ID
Index <- rowSums(expected_count_Total[, ..Indi1] >= 6) >= 0.2*length(Indi1) & rowSums(TPM_Total[, ..Indi1] >= 0.1) >= 0.2*length(Indi1)
ExpMat <- data.frame(TPM_Total[Index, ..Indi1])
row.names(ExpMat) <- TPM_Total[Index, gene_id]
stopifnot(identical(colnames(ExpMat), SampInfo_Total$ID))
colnames(ExpMat) <- SampInfo_Total$donor_id
TPM_Total <- ExpMat

fwrite(TPM_PolyA, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/PolyA_RNA_TPM.txt", row.names = TRUE, sep = "\t", quote = FALSE)
system("bgzip /mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/PolyA_RNA_TPM.txt")

fwrite(TPM_Total, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/Total_RNA_TPM.txt", row.names = TRUE, sep = "\t", quote = FALSE)
system("bgzip /mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/Total_RNA_TPM.txt")


# Transcript
library(data.table)
## PloyA RNA
SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/PolyA_RNA_raw_RSEM.RData")
Indi1 <- SampInfo_PolyA$ID
Index <- rowSums(expected_count_PolyA[, ..Indi1] >= 6) >= 0.2*length(Indi1) & rowSums(TPM_PolyA[, ..Indi1] >= 0.1) >= 0.2*length(Indi1)
ExpMat <- data.frame(TPM_PolyA[Index, ..Indi1])
row.names(ExpMat) <- TPM_PolyA[Index, transcript_id]
stopifnot(identical(colnames(ExpMat), SampInfo_PolyA$ID))
colnames(ExpMat) <- SampInfo_PolyA$donor_id
TPM_PolyA <- ExpMat

## Total RNA
SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt", sep = "\t")
load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/Total_RNA_raw_RSEM.RData")
Indi1 <- SampInfo_Total$ID
Index <- rowSums(expected_count_Total[, ..Indi1] >= 6) >= 0.2*length(Indi1) & rowSums(TPM_Total[, ..Indi1] >= 0.1) >= 0.2*length(Indi1)
ExpMat <- data.frame(TPM_Total[Index, ..Indi1])
sum(duplicated(TPM_Total[Index, transcript_id]))
row.names(ExpMat) <- TPM_Total[Index, transcript_id]
stopifnot(identical(colnames(ExpMat), SampInfo_Total$ID))
colnames(ExpMat) <- SampInfo_Total$donor_id
TPM_Total <- ExpMat

fwrite(TPM_PolyA, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/PolyA_RNA_TPM.txt", row.names = TRUE, sep = "\t", quote = FALSE)
system("bgzip /mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/PolyA_RNA_TPM.txt")

fwrite(TPM_Total, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/Total_RNA_TPM.txt", row.names = TRUE, sep = "\t", quote = FALSE)
system("bgzip /mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/transcript/Total_RNA_TPM.txt")



# Junction
library(data.table)
## PloyA RNA
SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/PolyA_RNA_Raw_SJs.RData")
SJ <- data.frame(SJ_tab[, -c(1:5)])
row.names(SJ) <- paste0(SJ_tab$seqnames, ":", SJ_tab$start, "-", SJ_tab$end, ":", SJ_tab$strand)
sum(duplicated(row.names(SJ)))
setkey(SampInfo_PolyA, ID)
colnames(SJ) <- SampInfo_PolyA[colnames(SJ), donor_id]

Index <- rowSums(SJ >= 1, na.rm = TRUE) >= 0.2*ncol(SJ) & rowSums(SJ, na.rm = TRUE) >= 10
SJ <- SJ[Index, ]
SJ[is.na(SJ)] <- 0
SJ_PolyA <- SJ


## Total RNA
SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt")
load(file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/Total_RNA_Raw_SJs.RData")
SJ <- data.frame(SJ_tab[, -c(1:5)])
row.names(SJ) <- paste0(SJ_tab$seqnames, ":", SJ_tab$start, "-", SJ_tab$end, ":", SJ_tab$strand)
sum(duplicated(row.names(SJ)))
setkey(SampInfo_Total, ID)
colnames(SJ) <- SampInfo_Total[colnames(SJ), donor_id]

Index <- rowSums(SJ >= 1, na.rm = TRUE) >= 0.2*ncol(SJ) & rowSums(SJ, na.rm = TRUE) >= 10
SJ <- SJ[Index, ]
SJ[is.na(SJ)] <- 0
SJ_Total <- SJ

fwrite(SJ_PolyA, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/PolyA_RNA_SJ.txt", row.names = TRUE, sep = "\t", quote = FALSE)
system("bgzip /mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/PolyA_RNA_SJ.txt")

fwrite(SJ_Total, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/Total_RNA_SJ.txt", row.names = TRUE, sep = "\t", quote = FALSE)
system("bgzip /mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/Total_RNA_SJ.txt")
