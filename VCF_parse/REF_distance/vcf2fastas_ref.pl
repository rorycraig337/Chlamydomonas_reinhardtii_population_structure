#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#takes a haploid VCF and will output a single FASTA file per individual
#usage: perl vcf2fastas.pl --vcf in.vcf

my $vcf;

GetOptions(
	'vcf=s' => \$vcf,
) or die;

open (IN, "$vcf") or die;

my @strains;

push @strains, "ref";

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^#CHROM/) {
	    my @splits = split(/\t/, $line);
	    my $split_count = 0;
	    foreach my $split (@splits) {
	        if ($split_count > 8) {
	        	push @strains, $split;
	        }
		$split_count++;
	    }
	}
	if ($line !~ /^#/) {
		last;
	}
}

close IN;

my $strain_number = 8;

foreach my $strain (@strains) {
	open (VCF, "$vcf") or die;
	open (OUT, ">$strain.fa") or die;
	print OUT ">$strain\n";
	my $line_count = 0;
	while (my $vcf_line = <VCF>) {
		chomp $vcf_line;
		if ($vcf_line !~ /^#/) {
			$line_count++;
			my @cols = split(/\t/, $vcf_line);
			my $ref = $cols[3];
			my $alt = $cols[4];
			if ($strain_number == 8) { 
				print OUT "$ref";
			}
			else {
				my @calls = split(/:/, $cols[$strain_number]);
				if ($calls[0] eq "0") {
					print OUT "$ref";
				}
				elsif ($calls[0] eq "1") {
					print OUT "$alt";
				}
				elsif ($calls[0] eq "null") {
					print OUT "-";
				}
				else {
					die;
				}
			}
			if ($line_count == 60) {
				print OUT "\n";
				$line_count = 0;
			}
		}
	}
	$strain_number++;
	unless ($line_count == 0) {
		print OUT "\n";
	}
	close VCF;
	close OUT;
}

exit;
