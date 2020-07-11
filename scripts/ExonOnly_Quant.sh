# totalRNA, exon,  strand

for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'total RNA'|awk -F ',' '{print $1}')
do
  htseq-count -m union -f bam -t exon -s reverse -r pos -i gene_id --nonunique all /mnt/data2/cell_line/EGAD00001002671_bam/${i}.Aligned.sortedByCoord.out.bam  /mnt/raid61/Personal_data/chenli/REF/genecode37.v30/gencode.v30lift37.annotation.gtf > ${i}".Exons.counts.total.txt"
done



# polyARNA, exon, non-strand

for i in $(cat /mnt/data1/chenli/cell_T/human_Tcell_2016_sampleinfo.csv|grep 'polyA RNA'|awk -F ',' '{print $1}'|awk '{if(NR >2 ){print}}')
do
  htseq-count -m union -f bam -t exon -s no -r pos -i gene_id --nonunique all /mnt/data2/cell_line/EGAD00001002671_bam/${i}.Aligned.sortedByCoord.out.bam  /mnt/raid61/Personal_data/chenli/REF/genecode37.v30/gencode.v30lift37.annotation.gtf > ${i}".Exons.counts.polyA.txt"
done


