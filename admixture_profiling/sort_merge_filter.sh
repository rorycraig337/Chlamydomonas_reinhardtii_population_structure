#!/usr/bin/env bash

for i in *introgressed.bed ; do
	sort -k1,1 -k2n,2n $i > s_$i
	rm $i
	bedtools merge -i s_$i > m_$i
	rm s_$i
	perl filter_bed_by_size.pl --bed m_$i --min_size 36000 --out $i
	rm m_$i
done
