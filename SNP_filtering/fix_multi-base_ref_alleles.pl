#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#if there is a deltion in a strain that has been removed from this VCF, the original REF allele is still >1 base in length
#script replaces string of bases with first base, which corresponds to the base at that position
#perl: fix_multi-base_ref_alleles.pl --vcf in.vcf --out out.vcf

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
	if ($line !~ /^#/) {
		my @cols = split(/\t/, $line);
		if (length $cols[3] > 1) {
			my $base = substr($cols[3], 0, 1);
			$cols[3] = $base;
			print OUT join ("\t", @cols);
			print OUT "\n";
			my $strain = 9;
			until ($strain == 34) {
				my @calls = split(/:/, $cols[$strain]);
				$strain++;
				if ($calls[0] eq "1") {
					die "$line\n";
				}
			}
		}
		else {
			print OUT "$line\n";
		}
	}
	else {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;
