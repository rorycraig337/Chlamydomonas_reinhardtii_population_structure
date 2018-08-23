#!/usr/bin/env bash

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R $PATH/chlamy.5.3.w_organelles_mtMinus.fasta -V $PATH/filter_snps.vcf --exclude_sample_file exclude_JPN.txt --removeUnusedAlternates -o JPN.snps.vcf 
