#!/usr/bin/env bash

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /localdisk/data/chlamydomonas/popgen/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/genome/chlamy.5.3.w_organelles_mtMinus.fasta -V /localdisk/home/s0920593/projects/biogeography/github/all_samples.vcf.gz --exclude_sample_file exclude.txt --removeUnusedAlternates -o species_wide.vcf

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /localdisk/data/chlamydomonas/popgen/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/genome/chlamy.5.3.w_organelles_mtMinus.fasta -V species_wide.vcf -selectType SNP -o species_wide.raw_snps.vcf

java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /localdisk/data/chlamydomonas/popgen/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/genome/chlamy.5.3.w_organelles_mtMinus.fasta -V species_wide.vcf -selectType INDEL -o species_wide.raw_indels.vcf
