#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to filter CDS annotation table for codons with more than one SNP (confusing designation of degeneracy)
#takes output from CDS_overlap_filter.pl, with SNP column modified by alter_allele_counts.pl 
#usage: perl multi-allelic_codon_filter.pl --table in.tsv --out out.tsv

my $table;
my $out;

GetOptions(
        'table=s' => \$table,
        'out=s' => \$out,
) or die "missing input\n";

open (IN, "$table") or die;
open (OUT, ">$out") or die;

my $line1;
my $line2;
my $line3;

my $pass;
my $snp_count;

my $codon_count = 0;
my $multi_count = 0;

while (my $line = <IN>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		my @strands = split(/:/, $cols[21]); #strand (+,-) information
		my @frames = split(/:/, $cols[22]); #frame (0,1,2) information
		if ($strands[0] eq "+") { #gene orientation
			if ($cols[13] == 1) { #coding sequence
				if ($frames[0] == 0) { #1st position
					$line1 = "$line";
					$snp_count = 0;
					$pass = 0;
					$pass++;
					my $bases = pop @cols; #assumes alleles are last column in line
      					my @alleles = split(/:/, $bases);
      					my $allele_count = 0;
      					foreach my $allele (@alleles) {
      						if ($allele > 0) {
      							$allele_count ++;
      						}
      					}
      					if ($allele_count > 1) {
      						$snp_count ++; #site is a snp
      					}
				}
				if ($frames[0] == 1) { #2nd position
					$line2 = "$line";
					$pass++;
					my $bases = pop @cols; #assumes alleles are last column in line
      					my @alleles = split(/:/, $bases);
      					my $allele_count = 0;
      					foreach my $allele (@alleles) {
      						if ($allele > 0) {
      							$allele_count ++;
      						}
      					}
      					if ($allele_count > 1) {
      						$snp_count ++; #site is a snp
     					}
				}
				if ($frames[0] == 2) { #3rd position
					$line3 = "$line";
					$pass++;
					my $bases = pop @cols; #assumes alleles are last column in line
      					my @alleles = split(/:/, $bases);
      					my $allele_count = 0;
      					foreach my $allele (@alleles) {
      						if ($allele > 0) {
      							$allele_count ++;
      						}
      					}
      					if ($allele_count > 1) {
      						$snp_count ++; #site is a snp
      					}
				}
				if ($pass == 3) { #contiguous condon (strand +, frame 0 -> 1 -> 2)
					$codon_count++;
					unless ($snp_count > 1) { #filters multi-snp codons
						print OUT "$line1\n$line2\n$line3\n";
					}
					else {
						$multi_count++;
					}
				}
			}
		}
		if ($strands[0] eq "-") { #gene orientation is reverse strand, reverse frame order
			if ($cols[13] == 1) { #coding sequence
				if ($frames[0] == 2) { #1st position
					$line1 = "$line";
					$snp_count = 0;
					$pass = 0;
					$pass++;
					my $bases = pop @cols; #assumes alleles are last column in line
      					my @alleles = split(/:/, $bases);
      					my $allele_count = 0;
      					foreach my $allele (@alleles) {
      						if ($allele > 0) {
      							$allele_count ++;
      						}
      					}
      					if ($allele_count > 1) {
      						$snp_count ++; #site is a snp
      					}
				}
				if ($frames[0] == 1) { #2nd position
					$pass++;
					$line2 = "$line";
					my $bases = pop @cols; #assumes alleles are last column in line
      					my @alleles = split(/:/, $bases);
      					my $allele_count = 0;
      					foreach my $allele (@alleles) {
      						if ($allele > 0) {
      							$allele_count ++;
      						}
      					}
      					if ($allele_count > 1) {
      						$snp_count ++;
      					}
				}
				if ($frames[0] == 0) { #3rd position
					$pass++;
					$line3 = "$line";
					my $bases = pop @cols; #assumes alleles are last column in line
      					my @alleles = split(/:/, $bases);
      					my $allele_count = 0;
      					foreach my $allele (@alleles) {
      						if ($allele > 0) {
      							$allele_count ++;
      						}
      					}
      					if ($allele_count > 1) {
      						$snp_count ++;
      					}
				}
				if ($pass == 3) { #contiguos codon (-, 2 -> 1 -> 0)
					$codon_count++;
					unless ($snp_count > 1) { #filters multi-snp codons
						print OUT "$line1\n$line2\n$line3\n";
					}
					else {
						$multi_count++;
					}
				}
			}
		}
	}
}

print "codons are $codon_count, multi-SNP are $multi_count\n";

close IN;
close OUT;

exit;
