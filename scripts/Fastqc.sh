## run fastqc
cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $1}'|xargs -P 3 -I {} fastqc -o /mnt/data1/chenli/cell_T/FastQC_result -f fastq /mnt/data1/liyupeng/data/project/EGAD00001002671_fq/{}.pair1.fastq.gz

cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $1}'|xargs -P 3 -I {} fastqc -o /mnt/data1/chenli/cell_T/FastQC_result -f fastq /mnt/data1/liyupeng/data/project/EGAD00001002671_fq/{}.pair2.fastq.gz

cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}'|xargs -P 3 -I {} fastqc -o /mnt/data1/chenli/cell_T/FastQC_result -f fastq /mnt/data1/liyupeng/data/project/EGAD00001002671_fq/{}.pair1.fastq.gz

cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}'|xargs -P 3 -I {} fastqc -o /mnt/data1/chenli/cell_T/FastQC_result -f fastq /mnt/data1/liyupeng/data/project/EGAD00001002671_fq/{}.pair2.fastq.gz



## run MultiQC
multiqc ./totalRNA/*_fastqc.zip -o ./multiQC -n multiqc_report_totalRNA  ##totalRNA
multiqc ./polyARNA/*_fastqc.zip -o ./multiQC -n multiqc_report_polyARNA  ##polyARNA

## extract base data(fastqc) to *base.txt (change grep and pair1/2)
for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $1}')
do
sed -n "13,68p" ${i}.pair1_fastqc/fastqc_data.txt > ${i}.pair1_base.txt
done

## extract duplicate data(fastqc) to *duplicate.txt (change grep and pair1/2)
for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $2}')
do
less tmp/${i}.pair2_fastqc/fastqc_data.txt |grep -A16 '#Duplication Level' >${i}.pair2_duplicate.txt
done


## extract GC%
# tmp <- fread('/mnt/data1/chenli/cell_T/FastQC_result/multiQC/multiqc_report_polyARNA_data/multiqc_fastqc.txt')
# tmp <- tmp[,.(GC = mean(`%GC`)),by=donor]
# write.table(tmp, file = "/mnt/data1/chenli/cell_T/FastQC_result/multiQC/multiqc_report_polyARNA_data/polyA_GC.txt", row.names = F, col.names=F,quote = F, sep="\t")

## calculate Q30
for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}')
do
python /mnt/data1/chenli/soft/q30/q30.py /mnt/raid61/Cell_2016_chen/EGAD00001002671/EGAD00001002671_fq_decrypted/${i}_*.pair2.fastq.gz > ${i}.Q30.pair2.txt
done

for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}')
do
python /mnt/data1/chenli/soft/q30/q30.py /mnt/raid61/Cell_2016_chen/EGAD00001002671/EGAD00001002671_fq_decrypted/${i}_*.pair1.fastq.gz > ${i}.Q30.pair1.txt
done


## extract Q30
for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}')
do
cat ${i}.Q30.pair1.txt|sed -n '5p'|awk -F '[,)]' '{print $2}' >> summary.p1.txt
done


for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}')
do
cat ${i}.Q30.pair2.txt|sed -n '5p'|awk -F '[,)]' '{print $2}' >> summary.p2.txt
done



