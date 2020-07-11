#### exon-only-use rlog

# sample infomation
SampInfo_PolyA <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/PolyA_RNA_sampleInfo.txt")
SampInfo_Total <- fread("/mnt/raid61/Personal_data/tangchao/ScientificData/data/SampleInfo/Total_RNA_sampleInfo.txt")
SampInfo <- rbind(SampInfo_PolyA,SampInfo_Total)
SampInfo <- tibble::column_to_rownames(setDF(SampInfo), var = 'ID')


# loading data
exonOnly_polyA <- read.table("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/Gene_Count_Only_Exon_PolyA.txt")
exonOnly_total <- read.table("/mnt/raid61/Personal_data/tangchao/ScientificData/data/Phenotype/gene/Gene_Count_Only_Exon_Total.txt")

exonOnly_cnt <- cbind(exonOnly_polyA,exonOnly_total)
exonOnly_cnt <- exonOnly_cnt[,rownames(SampInfo)]
exonOnly_cnt <- exonOnly_cnt[-((nrow(exonOnly_cnt)-4):nrow(exonOnly_cnt)),]

dds_EO <- DESeqDataSetFromMatrix(countData = exonOnly_cnt,
                                 colData = SampInfo,
                                 design = ~ 1)

dds_EO <- dds_EO[apply(assay(dds_EO),1,function(x){ sum(x >10)>16 }),]
rld_EO <- rlog(dds_EO, blind = FALSE)  # rlog



## combat
library(sva)
modcombat = model.matrix(~1, data = SampInfo)
batch = rep(c(1,2),each=40)
combat_rld_EO = ComBat(dat=assay(rld_EO), batch = batch, mod= modcombat)


##cluster
colnames(combat_rld_EO) <- SampInfo$Name
colnames(combat_rld_EO) <- c(gsub("_mRNA"," (PolyA-selected)",colnames(combat_rld_EO)[1:40]),
                             gsub("_RNA"," (rRNA-depleted)",colnames(combat_rld_EO)[41:80]))

pdf("EO_rlog_combat_dendrogram_rec2.pdf",width=13,height=12)
fviz_dend(hcut(t(combat_rld_EO), k = 4, stand = TRUE), cex=0.7,
          horiz=TRUE, color_labels_by_k=FALSE,label_cols="black",k_colors="black")
dev.off()