# Fig2a. Per Base Mean Sequence Quality

Base_Median = function(mypath,pattern){
  filenames = list.files(path = mypath, pattern = pattern, full.names = TRUE)
  datalist = lapply(filenames,function(x){
    tmp <- fread(x)
    colnames(tmp)<-c("Base","Mean","Median","Lower_Quartile","Upper_Quartile","Percentile_10th","Percentile_90th")
    tmp <-tmp[,c("Base","Median")]
    tmp$Base <- fct_inorder(factor(tmp$Base))  ## keep orders
    setnames(tmp,"Median",strsplit(tail(strsplit(x, "/")[[1]], n = 1),"[_]")[[1]][1])
    return(tmp)
  })
  Reduce(function(x, y) {merge(x, y, all = T, by = "Base")}, datalist)
}
system.time(Total_baseMedian<-Base_Median(mypath='/mnt/data1/chenli/cell_T/FastQC_result/totalRNA/',pattern='base.txt'))
system.time(PolyA_baseMedian<-Base_Median(mypath='/mnt/data1/chenli/cell_T/FastQC_result/polyARNA/',pattern='base.txt'))


## total
Total_baseMedian<-melt(Total_baseMedian)


## polyA
PolyA_baseMedian<-melt(PolyA_baseMedian)


setDF(Total_baseMedian)
setDF(PolyA_baseMedian)

all_base <- rbind(Total_baseMedian,PolyA_baseMedian)
all_base$library <- rep(c("totalRNA","polyARNA"),each=4400)



## plot
p1 <- ggplot(all_base, aes(x=Base, y=value, color=library)) +
  stat_boxplot(geom = "errorbar",width=0.8)+
  geom_boxplot(width=0.8,outlier.shape=NA)+
  labs(title="Per Base Mean Sequence Quality", x="Position (bp)", y = "Mean Sequence Quality (Phred Score)")+
  theme_classic() + 
  scale_color_manual(values=c("#00AFBB", "#E7B800"))+
  #scale_y_continuous(breaks=c(0,10,20,30,40),labels=c("0","10","20","30","40"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.text=element_text(size=10),
        legend.title= element_blank())+
  expand_limits(y=c(0,41))+scale_y_continuous(expand = c(0,0)) 





# Fig2b.genebody coverage

file_list <- list.files("/mnt/data1/chenli/cell_T/RSeQC_result/gene_body", pattern = "geneBodyCoverage.txt")
file_list <- paste0("/mnt/data1/chenli/cell_T/RSeQC_result/gene_body/",file_list)
data <- lapply(file_list,function(x){
  tab <- fread(x,sep = '\t', header = F)
  tab <- tab[2,]
})
geneBodyCoverage <- do.call(rbind,data)
geneBodyCoverage2<- apply(geneBodyCoverage[,-1],1,function(x){
  (x - min(x)) / (max(x) - min(x))   ## txt files  (i - min(dat)) / (max(dat) - min(dat))
})
colnames(geneBodyCoverage2)<- do.call(rbind, strsplit(geneBodyCoverage$V1,split="[.]"))[,1] 
geneBodyCoverage2 <- cbind(geneBodyCoverage2,percent=1:100)
geneBodyCoverage2 <- melt(geneBodyCoverage2,id.var="percent")
geneBodyCoverage2 <- geneBodyCoverage2[-(8001:8100),]
geneBodyCoverage2$Var1<-rep(1:100,80)
geneBodyCoverage2$library <- rep(c(rep("polyARNA",100),rep("totalRNA",100)),40)


p2 <- ggplot(geneBodyCoverage2,aes(x=Var1,y=value,group=Var2))+geom_line(aes(colour=library))+
          theme(panel.background = element_blank(), 
                axis.title = element_text(size = 11), 
                axis.text = element_text(size = 10), 
                axis.line = element_line(), 
                title = element_text(size = 15, hjust = .5, vjust = .5), 
                plot.title = element_text(hjust = 0.5), 
                plot.subtitle = element_text(hjust = 0.5),
                legend.position="none")+
          scale_colour_manual(values=c("#00AFBB", "#E7B800"))+
          labs(x="Gene body percentile (5'->3')",y="Coverage",title="Gene body coverage graph")






#Fig2c.The percentages of counts in various gene regions

sampleinfo <-read.csv('/mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv')

RSeqc_summary <- function(column=c(1,2),pattern=".UTR.txt",by="Group"){
  filenames2 <- list.files("/mnt/data1/chenli/cell_T/RSeQC_result/read_distribution/other", pattern = pattern)
  filenames2 <- paste0("/mnt/data1/chenli/cell_T/RSeQC_result/read_distribution/other/",filenames2)
  RSeqc_file <-lapply(filenames2,function(x){
    tmp <-fread(x,select=column)
    colnames(tmp)[2] <- strsplit(tail(strsplit(x, "/")[[1]], n = 1),"[.]")[[1]][1]
    return(tmp)
  })
  Reduce(function(x, y) {merge(x, y, all = T, by = by)}, RSeqc_file) 
}

RSeqc_file_totalbases <- RSeqc_summary(column=c(1,2))
RSeqc_file_Tag_count <- RSeqc_summary(column=c(1,3))
RSeqc_file_TagsperKb <- RSeqc_summary(column=c(1,4))

colnames(RSeqc_file_totalbases)[-1] <-sampleinfo$Name
colnames(RSeqc_file_Tag_count)[-1]  <-sampleinfo$Name
colnames(RSeqc_file_TagsperKb)[-1]  <-sampleinfo$Name


Tag_count <- melt(RSeqc_file_Tag_count[c(1:5,8),])
Tag_count2 <- Tag_count[,.(percentage = value*100/sum(value)),by=variable]
Tag_count$percentage <-Tag_count2$percentage
Tag_count$Library <-rep(rep(c("polyARNA","totalRNA"),each=6),40)
df <- Tag_count[,.(mean = mean(percentage),sd=sd(percentage)),by=.(Library,Group)]


p3 <- ggplot(df, aes(x=Group, y=mean, fill=Library)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  labs(y = "Frequency (%)")+
  theme_classic() +
  scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
  theme(axis.title = element_text(size = 16),
        axis.title.x=element_blank(),
        axis.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text=element_text(size=10),
        legend.title= element_blank())