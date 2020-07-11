library(GenomicFeatures)
library(data.table)
gtfFile <- file.path("/mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/gencode.v30lift37.annotation.gtf")
txdb <- makeTxDbFromGFF(file=gtfFile,
                        format="gtf",
                        dataSource="ensembl",
                        organism="Homo sapiens")
isActiveSeq(txdb)
ebg <- exonsBy(txdb, by="gene")
ebg

## Total RNA
SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt", sep = "\t")
bamFile <- paste("/mnt/raid62/Chen_Cell_2016/EGAD00001002671_bam/", SampInfo_Total$ID , ".Aligned.sortedByCoord.out.bam", sep = "")
stopifnot(all(file.exists(bamFile)))

library(GenomicAlignments)
library(BiocParallel)

flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE,
                    isUnmappedQuery = FALSE,
                    isDuplicate = NA)
sbp <- ScanBamParam(flag=flag, mapqFilter = 255)

bamLst <- BamFileList(bamFile, yieldSize=2000000)
options(srapply_fapply = "parallel", mc.cores=40)
se_Total <- summarizeOverlaps(features = ebg, 
                              reads = bamLst, 
                              mode = "IntersectionStrict",
                              ignore.strand = FALSE, 
                              inter.feature = FALSE, 
                              singleEnd = FALSE,
                              fragments = FALSE, 
                              strandMode = 2, 
                              param = sbp, 
                              preprocess.reads = NULL)

## PolyA RNA

SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
bamFile <- paste("/mnt/raid62/Chen_Cell_2016/EGAD00001002671_bam/", SampInfo_PolyA$ID , ".Aligned.sortedByCoord.out.bam", sep = "")
stopifnot(all(file.exists(bamFile)))

bamLst <- BamFileList(bamFile, yieldSize=2000000)
options(srapply_fapply = "parallel", mc.cores=20)
se_PolyA <- summarizeOverlaps(features = ebg, 
                              reads = bamLst, 
                              mode = "IntersectionStrict",
                              ignore.strand = TRUE, 
                              inter.feature = FALSE, 
                              singleEnd = FALSE,
                              fragments = FALSE, 
                              param = sbp, 
                              preprocess.reads = NULL)

save(se_Total, se_PolyA, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/se_UnrmDuplicate.RData")
