## genebody coverage
cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'RNA'|awk -F ',' '{print $1}'|xargs -P 3 -I {} python /mnt/data1/chenli/soft/geneBody_coverage.py \
 -i /mnt/data2/cell_line/EGAD00001002671_bam/{}.Aligned.sortedByCoord.out.bam \
 -r /mnt/raid61/Personal_data/chenli/other/gencode.v30lift37.annotation.gtf.bed12 \
 -o {}

## reads distribution
for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $1}'|awk '{if(NR <21 ){print}}')
do
read_distribution.py -i /mnt/data2/cell_line/EGAD00001002671_bam/${i}.Aligned.sortedByCoord.out.bam -r /mnt/data1/chenli/mydata/gencode_37/gencode.v30lift37.annotation.gtf.bed12 > ${i}_result.txt
done

for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $1}'|awk '{if(NR >20 ){print}}')
do
read_distribution.py -i /mnt/data2/cell_line/EGAD00001002671_bam/${i}.Aligned.sortedByCoord.out.bam -r /mnt/data1/chenli/mydata/gencode_37/gencode.v30lift37.annotation.gtf.bed12 > ${i}_result.txt
done

for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}'|awk '{if(NR >20 ){print}}')
do
read_distribution.py -i /mnt/data2/cell_line/EGAD00001002671_bam/${i}.Aligned.sortedByCoord.out.bam -r /mnt/data1/chenli/mydata/gencode_37/gencode.v30lift37.annotation.gtf.bed12 > ${i}_result.txt
done

for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}'|awk '{if(NR <21 ){print}}')
do
read_distribution.py -i /mnt/data2/cell_line/EGAD00001002671_bam/${i}.Aligned.sortedByCoord.out.bam -r /mnt/data1/chenli/mydata/gencode_37/gencode.v30lift37.annotation.gtf.bed12 > ${i}_result.txt
done

## extract various gene regions
for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'RNA'|awk -F ',' '{print $1}')
do
sed -n "5,15p" /mnt/data1/chenli/cell_T/RSeQC_result/read_distribution/${i}_result.txt > ${i}.UTR.txt
done

## TIN (transcript integrity number)
cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $1}' | xargs -P 10 -I {} tin.py -i /mnt/data2/cell_line/EGAD00001002671_bam/{}.Aligned.sortedByCoord.out.bam -r /mnt/raid61/Personal_data/chenli/other/gencode.v30lift37.annotation.gtf.bed12 -s

cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}'|xargs -P 10 -I {} tin.py -i /mnt/data2/cell_line/EGAD00001002671_bam/{}.Aligned.sortedByCoord.out.bam -r /mnt/raid61/Personal_data/chenli/other/gencode.v30lift37.annotation.gtf.bed12


cat *.summary.txt|awk '!(NR%2){print}' >TIN_summary.txt

