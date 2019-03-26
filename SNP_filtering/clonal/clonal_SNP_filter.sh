#!/usr/bin/env bash

gunzip ../../VCF_parse/clonal/*vcf.gz

for i in ../../VCF_parse/clonal/*vcf ; do
	root=`basename $i .vcf`
	#filter indels
	perl ../indel_filter_and_rescue.pl --vcf $i --out $root.f1.vcf --flank 5 --suffix indel_mask.bed --mode rescue
	bash ../sort_and_merge.indels.sh
	perl ../../differentiation_divergence_analyses/mask_vcf_regions.pl --samples CC1009,CC1010,CC1373,CC1952,CC2342,CC2343,CC2344,CC2931,CC2932,CC2935,CC2936,CC2937,CC2938,CC3059,CC3060,CC3061,CC3062,CC3063,CC3064,CC3065,CC3068,CC3069,CC3071,CC3073,CC3075,CC3076,CC3079,CC3082,CC3083,CC3084,CC3086,CC3268,GB117,GB119,GB123,GB13,GB138,GB141,GB57,GB66,NIES2463,NIES2464 --vcf $root.f1.vcf --suffix .indel_mask.bed --out $root.f2.vcf
	perl ../SNP_validator.pl --vcf $root.f2.vcf --out $root.f3.vcf
	#filter on max nuclear depth and quality
	perl ../failed_genotype_filter_maxDP.pl --vcf $root.f3.vcf --GQ 20 --DP_min 3 --DP_max ../../VCF_parse/raw_coverage/nuclear_cov.txt --DP_max_int 4 --out $root.f4.vcf
	perl ../SNP_validator.pl --vcf $root.f4.vcf --out $root.f5.vcf
	#extract SNPs
	perl ../extract_SNPs.pl --vcf $root.f5.vcf --out snps.$root.vcf
	rm *.f*.vcf
	rm *indel_mask.bed
	line=$(cat snps.$root.vcf | wc -l)
	hash=$(grep -c '#' snps.$root.vcf)
	total=$(($line - $hash))
	echo "$root SNP count: $total"
done

gzip ../../VCF_parse/clonal/*vcf
