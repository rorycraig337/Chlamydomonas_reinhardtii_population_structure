#!/usr/bin/env bash

for i in *IBD.bed ; do
        sort -k1,1 -k2n,2n $i > s_$i
        rm $i
	bedtools merge -i s_$i > m_$i
        rm s_$i
        mv m_$i $i
        cat $i | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
