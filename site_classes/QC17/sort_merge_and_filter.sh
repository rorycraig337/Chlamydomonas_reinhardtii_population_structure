#!/usr/bin/env bash

for i in *bed ; do
	root=`basename $i .bed`
	sort -k1,1 -k2n,2n $i | bedtools merge -i stdin > s_$i
	grep -v "mtMinus" s_$i | grep -v "mtDNA" | grep -v "cpDNA" | bedtools subtract -a stdin -b ../../SNP_filtering/MT.bed | bedtools subtract -a stdin -b ../Ns.bed3 > QC17_$i
	bedtools getfasta -fi /localdisk/data/chlamydomonas/popgen/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/genome/chlamy.5.3.w_organelles_mtMinus.fasta -bed QC17_$i > QC17_$root.fa
	perl get_GC.pl --fasta QC17_$root.fa
	rm s_$i
        echo "$i:"
        cat QC17_$i | awk 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
