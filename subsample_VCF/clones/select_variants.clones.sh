#!/usr/bin/env bash

for i in clone*txt ; do
	root=`basename $i .txt`
	java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R $PATH/chlamy.5.3.w_organelles_mtMinus.fasta -V $PATH/all_samples.vcf.gz --sample_file $i --removeUnusedAlternates -selectType SNP -o $root.raw_snps.vcf
	java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R $PATH/chlamy.5.3.w_organelles_mtMinus.fasta -V $PATH/all_samples.vcf.gz --sample_file $i --removeUnusedAlternates -selectType INDEL -o $root.raw_indels.vcf
done
