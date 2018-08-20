#!/usr/bin/env bash

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R $PATH/chlamy.5.3.w_organelles_mtMinus.fasta -V $PATH/all_samples.vcf.gz --exclude_sample_file exclude.txt --removeUnusedAlternates -o species_wide.vcf

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R $PATH/chlamy.5.3.w_organelles_mtMinus.fasta -V species_wide.vcf -selectType SNP -o species_wide.raw_snps.vcf

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R $PATH/chlamy.5.3.w_organelles_mtMinus.fasta -V species_wide.vcf -selectType INDEL -o species_wide.raw_indels.vcf
