#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to convert VCF file of SNPs to input file for hmmIBD
#usage: perl vcf_to_hmmIBD.pl --vcf snps.vcf --out snps.hmmIBD.txt

my $vcf;
my $out;

GetOptions(
        'vcf=s' => \$vcf,
        'out=s' => \$out,
) or die "missing input\n";

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^#CHROM/) {
		my @head_cols = split(/\t/, $line);
		my $head_col_count = 0;
		print OUT "chrom\tpos";
		foreach my $head_col (@head_cols) {
			if ($head_col_count > 8) {
				print OUT "\t$head_col"; #print header of isolate IDS	
			}
			$head_col_count++;
		}
		print OUT "\n";
	}
	unless ($line =~ /^#/) {
		if ($line =~ /scaff/) {
			next; #ignore unassembled scaffolds
		}
		my @cols = split(/\t/, $line);
		my $chrom = substr $cols[0], 11; #reduce chromosome to number only
		my $col_count = 0;
		print OUT "$chrom\t$cols[1]";
		foreach my $col (@cols) {
			if ($col_count > 8) {
				my @calls = split(/:/, $col);
				if ($calls[0] eq "0") {
					print OUT "\t0";
				}
				elsif ($calls[0] eq "1") {
					print OUT "\t1";
				}
				else {
					print OUT "\t-1"; #missing data
				}
			}
			$col_count++;
		}
		print OUT "\n";	
	}
}

close IN;
close OUT;

exit;
