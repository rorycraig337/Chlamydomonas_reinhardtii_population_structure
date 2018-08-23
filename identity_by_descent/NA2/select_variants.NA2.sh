#!/usr/bin/env bash

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R $PATH/chlamy.5.3.w_organelles_mtMinus.fasta -V $PATH/filter_snps.vcf --exclude_sample_file exclude_NA2.txt --removeUnusedAlternates -o NA2.snps.vcf 
