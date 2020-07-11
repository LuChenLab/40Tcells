library(SCDSS)
library(data.table)

## PolyA ----

mypath <- file.path("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/SJ/PolyA")

SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
sampinfo <- setDF(SampInfo_PolyA, rownames = SampInfo_PolyA[[1]])
gtfFile <- file.path("/mnt/raid61/Personal_data/tangchao/Document/gencode/human/GRCh37/gencode.v30lift37.annotation.gtf")

mysetFromSJ <- SCASDataSetFromSJFiles(path = mypath,
                                      pattern = ".SJ.out.tab",
                                      SampleInfo = sampinfo,
                                      design = ~ 1,
                                      minSJ = 2,
                                      SSMinSJ = 10,
                                      minSJs = 100,
                                      minSamps = 2,
                                      uniqueMapOnly = TRUE,
                                      gtfFile = gtfFile,
                                      NT = 4,
                                      project.name = "PolyA")
mysetFromSJ

library(SummarizedExperiment)
counts <- assay(mysetFromSJ, 1)
row.names(counts) <- as.character(rowRanges(mysetFromSJ))
start.psi <- assay(mysetFromSJ, 2)
row.names(start.psi) <- paste0(as.character(rowRanges(mysetFromSJ)), "|ss")
end.psi <- assay(mysetFromSJ, 3)
row.names(end.psi) <- paste0(as.character(rowRanges(mysetFromSJ)), "|se")
sjInfo <- rowRanges(mysetFromSJ)

colnames(counts) <- sampinfo[colnames(counts), ]$donor_id
colnames(start.psi) <- sampinfo[colnames(start.psi), ]$donor_id
colnames(end.psi) <- sampinfo[colnames(end.psi), ]$donor_id

PSI <- rbind(start.psi, end.psi)
PSI <- PSI[which(rowSds(PSI, na.rm = T) > 0), ]
PSI <- data.frame(PSI)

saveRDS(PSI, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/PolyA_intron-centric_PSI.Rds")

PolyA_PSI <- copy(PSI)


## Total ----

mypath <- file.path("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/SJ/Total/")

SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt")
sampinfo <- setDF(SampInfo_Total, rownames = SampInfo_Total[[1]])

mysetFromSJ <- SCASDataSetFromSJFiles(path = mypath,
                                      pattern = ".SJ.out.tab",
                                      SampleInfo = sampinfo,
                                      design = ~ 1,
                                      minSJ = 2,
                                      SSMinSJ = 10,
                                      minSJs = 100,
                                      minSamps = 2,
                                      uniqueMapOnly = TRUE,
                                      gtfFile = gtfFile,
                                      NT = 4,
                                      project.name = "Total")
mysetFromSJ

counts <- assay(mysetFromSJ, 1)
row.names(counts) <- as.character(rowRanges(mysetFromSJ))
start.psi <- assay(mysetFromSJ, 2)
row.names(start.psi) <- paste0(as.character(rowRanges(mysetFromSJ)), "|ss")
end.psi <- assay(mysetFromSJ, 3)
row.names(end.psi) <- paste0(as.character(rowRanges(mysetFromSJ)), "|se")
sjInfo <- rowRanges(mysetFromSJ)

colnames(counts) <- sampinfo[colnames(counts), ]$donor_id
colnames(start.psi) <- sampinfo[colnames(start.psi), ]$donor_id
colnames(end.psi) <- sampinfo[colnames(end.psi), ]$donor_id

PSI <- rbind(start.psi, end.psi)
PSI <- PSI[which(rowSds(PSI, na.rm = T) > 0), ]
PSI <- data.frame(PSI)

saveRDS(PSI, "/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/junction/Total_intron-centric_PSI.Rds")
