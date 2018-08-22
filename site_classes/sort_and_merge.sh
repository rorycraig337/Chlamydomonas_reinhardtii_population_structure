#!/usr/bin/env bash

for i in *bed ; do
	sort -k1,1 -k2n,2n $i | bedtools merge -i stdin > s_$i
	mv s_$i $i
done
