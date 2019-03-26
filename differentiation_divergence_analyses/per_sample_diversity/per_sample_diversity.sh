#extract isolates from VCF

bedtools intersect -a ../../SNP_filtering/species_wide/all_sites.min_isolates.nuclear.1lab.vcf -b ../../site_classes/4D.bed -header > 4D.vcf

for i in *_ex.txt ; do
	root=`basename $i _ex.txt`
	perl ../../SNP_filtering/vcf_isolate_filter.pl --in 4D.vcf --exclude $i --out $root.vcf
	perl ../../SNP_filtering/SNP_validator.pl --vcf $root.vcf --out $root.4D.vcf
done

perl ../chlamy_popgen.pl ../../site_classes/4D.bed NC.4D.vcf NC.4D --input_format vcf --bed_format BED6 --mode overall --output pi,TajD --min_copies 2 --max_copies 3
perl ../chlamy_popgen.pl ../../site_classes/4D.bed QC.4D.vcf QC.4D --input_format vcf --bed_format BED6 --mode overall --output pi,TajD --min_copies 2 --max_copies 24
perl ../chlamy_popgen.pl ../../site_classes/4D.bed MA.4D.vcf MA.4D --input_format vcf --bed_format BED6 --mode overall --output pi --min_copies 2 --max_copies 2
perl ../chlamy_popgen.pl ../../site_classes/4D.bed Far93.4D.vcf Far93.4D --input_format vcf --bed_format BED6 --mode overall --output pi,TajD --min_copies 2 --max_copies 14
perl ../chlamy_popgen.pl ../../site_classes/4D.bed Mac.4D.vcf Mac.4D --input_format vcf --bed_format BED6 --mode overall --output pi,TajD --min_copies 2 --max_copies 3
perl ../chlamy_popgen.pl ../../site_classes/4D.bed Far16.4D.vcf Far16.4D --input_format vcf --bed_format BED6 --mode overall --output pi,TajD --min_copies 2 --max_copies 7

