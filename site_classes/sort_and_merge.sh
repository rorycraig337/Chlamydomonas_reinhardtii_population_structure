#!/usr/bin/env bash

for i in *.bed ; do
	sort -k1,1 -k2n,2n $i | bedtools merge -i stdin > s_$i
	bedtools subtract -a s_$i -b Ns.bed3 > $i
	echo "$i:"
	cat $i | awk 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'	
	rm s_$i
done
