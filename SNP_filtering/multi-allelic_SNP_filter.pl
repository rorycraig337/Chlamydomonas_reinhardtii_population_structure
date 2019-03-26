#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#script filters SNPs that are multi-allelic (comma in 5th column) or that have a deletion (* in 5th column)
#usage: perl multi-allelic_SNP_filter.pl --vcf in.vcf --out out.vcf

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
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		unless ( ($cols[4] =~ /\*/) or ($cols[4] =~ /,/) ) {
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
