library(GenomicFeatures)
library(data.table)
gtfFile <- file.path("/mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/gencode.v30lift37.annotation.gtf")
txdb <- makeTxDbFromGFF(file=gtfFile,
                        format="gtf",
                        dataSource="ensembl",
                        organism="Homo sapiens")
isActiveSeq(txdb)
# isActiveSeq(txdb)[-match(paste0("chr", c(1:22, "X", "Y")), names(isActiveSeq(txdb)))] <- FALSE
names(which(isActiveSeq(txdb)))

gene_gr <- genes(txdb, columns="gene_id")

ebg <- exonsBy(txdb, by="gene")
ebg


SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt", sep = "\t")
bamFile <- paste("/mnt/raid62/Chen_Cell_2016/EGAD00001002671_bam/", SampInfo_Total$ID , ".Aligned.sortedByCoord.out.bam", sep = "")
stopifnot(all(file.exists(bamFile)))

library(parallel)
library(data.table)
library(GenomicAlignments)

# ----
# strand(gene_gr) <- ifelse(strand(gene_gr) == "+", "-", "+")

options(srapply_fapply="parallel", mc.cores=12)

# galp0 <- readGAlignmentPairs(bamFile, use.names = TRUE, param = sbp, strandMode = 2)

flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE,
                    isUnmappedQuery = FALSE,
                    isDuplicate = FALSE)
sbp <- ScanBamParam(flag=flag, mapqFilter = 255)

bamLst <- BamFileList(bamFile[1], yieldSize = 1000000)

SO1 <- summarizeOverlaps(features = gene_gr, 
                         reads = bamLst, 
                         mode = "IntersectionStrict", 
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE, 
                         fragments = TRUE, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO2 <- summarizeOverlaps(features = gene_gr, 
                         reads = bamLst, 
                         mode = "IntersectionStrict", 
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE, 
                         fragments = FALSE, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO3 <- summarizeOverlaps(features = gene_gr, 
                         reads = bamLst, 
                         mode = "IntersectionStrict", 
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE, 
                         fragments = FALSE, 
                         strandMode = 2, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO4 <- summarizeOverlaps(features = ebg, 
                         reads = bamLst, 
                         mode = "IntersectionStrict",
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE,
                         fragments = FALSE, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO5 <- summarizeOverlaps(features = ebg, 
                         reads = bamLst, 
                         mode = "IntersectionStrict",
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE,
                         fragments = FALSE, 
                         strandMode = 2, 
                         param = sbp, 
                         preprocess.reads = NULL)




SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
bamFile <- paste("/mnt/raid62/Chen_Cell_2016/EGAD00001002671_bam/", SampInfo_PolyA$ID , ".Aligned.sortedByCoord.out.bam", sep = "")
stopifnot(all(file.exists(bamFile)))
bamLst <- BamFileList(bamFile[1], yieldSize = 1000000)

SO6 <- summarizeOverlaps(features = gene_gr, 
                         reads = bamLst, 
                         mode = "IntersectionStrict", 
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE, 
                         fragments = TRUE, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO7 <- summarizeOverlaps(features = gene_gr, 
                         reads = bamLst, 
                         mode = "IntersectionStrict", 
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE, 
                         fragments = FALSE, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO8 <- summarizeOverlaps(features = gene_gr, 
                         reads = bamLst, 
                         mode = "IntersectionStrict", 
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE, 
                         fragments = FALSE, 
                         strandMode = 2, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO9 <- summarizeOverlaps(features = ebg, 
                         reads = bamLst, 
                         mode = "IntersectionStrict",
                         ignore.strand = FALSE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE,
                         fragments = FALSE, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO10 <- summarizeOverlaps(features = ebg, 
                          reads = bamLst, 
                          mode = "IntersectionStrict",
                          ignore.strand = FALSE, 
                          inter.feature = FALSE, 
                          singleEnd = FALSE,
                          fragments = FALSE, 
                          strandMode = 2, 
                          param = sbp, 
                          preprocess.reads = NULL)

SO11 <- summarizeOverlaps(features = gene_gr, 
                          reads = bamLst, 
                          mode = "IntersectionStrict", 
                          ignore.strand = TRUE, 
                          inter.feature = FALSE, 
                          singleEnd = FALSE, 
                          fragments = TRUE, 
                          param = sbp, 
                          preprocess.reads = NULL)

SO12 <- summarizeOverlaps(features = gene_gr, 
                         reads = bamLst, 
                         mode = "IntersectionStrict", 
                         ignore.strand = TRUE, 
                         inter.feature = FALSE, 
                         singleEnd = FALSE, 
                         fragments = FALSE, 
                         param = sbp, 
                         preprocess.reads = NULL)

SO13 <- summarizeOverlaps(features = ebg, 
                          reads = bamLst, 
                          mode = "IntersectionStrict",
                          ignore.strand = TRUE, 
                          inter.feature = FALSE, 
                          singleEnd = FALSE,
                          fragments = FALSE, 
                          param = sbp, 
                          preprocess.reads = NULL)

SO14 <- summarizeOverlaps(features = ebg, 
                          reads = bamLst, 
                          mode = "IntersectionStrict",
                          ignore.strand = TRUE, 
                          inter.feature = FALSE, 
                          singleEnd = FALSE,
                          fragments = TRUE, 
                          param = sbp, 
                          preprocess.reads = NULL)

save(SO1, SO2, SO3, SO4, SO5, SO6, SO7, SO8, SO9, SO10, SO11, SO12, SO13, SO14, file = "/mnt/raid61/Personal_data/tangchao/ScientificData/analysis/SO.test.RData")

# 5 与 13 最佳