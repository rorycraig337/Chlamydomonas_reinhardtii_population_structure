#!/usr/bin/env bash

for i in clone*txt ; do
	root=`basename $i .txt`
	java -jar /opt/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /localdisk/data/chlamydomonas/popgen/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/genome/chlamy.5.3.w_organelles_mtMinus.fasta -V /localdisk/home/s0920593/projects/biogeography/github/all_samples.vcf.gz --sample_file $i --removeUnusedAlternates -o $root.vcf
done
